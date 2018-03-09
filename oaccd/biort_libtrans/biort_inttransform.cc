/*
 * @BEGIN LICENSE
 *
 * oaccd by Psi4 Rolf H. Myhre April 2017, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "biort_inttransform.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libciomr/libciomr.h"

using namespace std;

namespace psi{ namespace oaccd {
BiortIntTransform::BiortIntTransform(std::shared_ptr<Wavefunction> wfn,
                                     SpaceVec spaces,
                                     SharedMatrix lCa,
                                     SharedMatrix rCa,
                                     SharedMatrix lCb,
                                     SharedMatrix rCb,
                                     TransformationType transformationType,
                                     OutputType outputType,
                                     MOOrdering moOrdering,
                                     FrozenOrbitals frozenOrbitals,
                                     bool init):
    lCa_(lCa),
    rCa_(rCa),
    lCb_(lCb),
    rCb_(rCb),
    IntegralTransform(wfn,
                      spaces,
                      transformationType,
                      outputType,
                      moOrdering,
                      frozenOrbitals,
                      init)
{
    // Implement set/get functions to customize any of this stuff.  Delayed initialization
    // is possible in case any of these variables need to be changed before setup.
    memory_ = Process::environment.get_memory();

    labels_  = wfn->molecule()->irrep_labels();
    nirreps_ = wfn->nirrep();
    nmo_     = wfn->nmo();
    nso_     = wfn->nso();
    sopi_    = wfn->nsopi();
    mopi_    = wfn->nmopi();
    clsdpi_  = wfn->doccpi();
    openpi_  = wfn->soccpi();
    frzcpi_  = wfn->frzcpi();
    frzvpi_  = wfn->frzvpi();
    nalphapi_ = wfn->nalphapi();
    nbetapi_ = wfn->nbetapi();
    frozen_core_energy_ = 0.0; 

    outfile->Printf("The Legend of Zelda \n");

    common_initialize();

    if(init) initialize();

}

/**
 * Sets up the DPD buffers and performs semicanonicalization, if necessary.
 */
void BiortIntTransform::initialize()
{
    print_         = Process::environment.options.get_int("PRINT");
    printTei_      = print_ > 5;
    useIWL_        = outputType_ == OutputType::IWLAndDPD || outputType_ == OutputType::IWLOnly;
    useDPD_        = outputType_ == OutputType::IWLAndDPD || outputType_ == OutputType::DPDOnly;
    iwlAAIntFile_  = transformationType_ == TransformationType::Restricted ? PSIF_MO_TEI : PSIF_MO_AA_TEI;
    iwlABIntFile_  = transformationType_ == TransformationType::Restricted ? PSIF_MO_TEI : PSIF_MO_AB_TEI;
    iwlBBIntFile_  = transformationType_ == TransformationType::Restricted ? PSIF_MO_TEI : PSIF_MO_BB_TEI;

    outfile->Printf("\n Heart of the Swarm \n\n");

    tpdm_buffer_ = 0;

    aQT_ = init_int_array(nmo_);
    if(transformationType_ == TransformationType::Restricted){
        reorder_qt(clsdpi_, openpi_, frzcpi_, frzvpi_, aQT_, mopi_, nirreps_);
        bQT_ = aQT_;
    }else{
        bQT_ = init_int_array(nmo_);
        reorder_qt_uhf(clsdpi_, openpi_, frzcpi_, frzvpi_, aQT_, bQT_, mopi_, nirreps_);
    }
    // Set up the correlated to Pitzer arrays.  These have to include the occupied core terms, because
    // the reference contributions are already folded into the TPDM.  However, they don't include frozen
    // virtuals
    aCorrToPitzer_ = init_int_array(nmo_);
    if(transformationType_ != TransformationType::Restricted){
        bCorrToPitzer_ = init_int_array(nmo_);
    }else{
        bCorrToPitzer_ = aCorrToPitzer_;
    }

    int nFzvFound = 0;
    int pitzerCount = 0;
    for(int h = 0; h < nirreps_; ++h){
        for (int p = 0; p < mopi_[h]; p++) {
            if (p < mopi_[h] - frzvpi_[h]) {
                // This is active, count it
                int q = aQT_[pitzerCount];
                aCorrToPitzer_[q] = pitzerCount - nFzvFound;
                if(transformationType_ != TransformationType::Restricted){
                    int q = bQT_[pitzerCount];
                    bCorrToPitzer_[q] = pitzerCount - nFzvFound;
                }
            }else{
                nFzvFound++;
            }
            pitzerCount++;
        }
    }

    if(print_ > 4){
        outfile->Printf( "\tThe Alpha Pitzer to QT mapping array:\n\t\t");
        for(int p = 0; p < nmo_; ++p)
            outfile->Printf( "%d ", aQT_[p]);
        outfile->Printf( "\n");
        outfile->Printf( "\tThe Beta Pitzer to QT mapping array:\n\t\t");
        for(int p = 0; p < nmo_; ++p)
            outfile->Printf( "%d ", bQT_[p]);
        outfile->Printf( "\n");
        outfile->Printf( "\tThe Alpha Correlated to Pitzer mapping array:\n\t\t");
        for(int p = 0; p < nmo_; ++p)
            outfile->Printf( "%d ", aCorrToPitzer_[p]);
        outfile->Printf( "\n");
        outfile->Printf( "\tThe Beta Correlated to Pitzer mapping array:\n\t\t");
        for(int p = 0; p < nmo_; ++p)
            outfile->Printf( "%d ", bCorrToPitzer_[p]);
        outfile->Printf( "\n");
    }

    process_spaces();

    // Set up the DPD library
    // TODO implement caching of files
    int numSpaces = spacesUsed_.size();
    int numIndexArrays = numSpaces * (numSpaces - 1) + 5 * numSpaces;
    cacheFiles_ = init_int_array(PSIO_MAXUNIT);
    cacheList_  = init_int_matrix(numIndexArrays, numIndexArrays);
    int currentActiveDPD = psi::dpd_default;
    dpd_init(myDPDNum_, nirreps_, memory_, 0, cacheFiles_, cacheList_, NULL, numSpaces, spaceArray_);

    // We have to redefine the MO coefficients for a UHF-like treatment
    if(transformationType_ == TransformationType::SemiCanonical){
        throw PSIEXCEPTION("Semicanonical is deprecated in Libtrans. Please pre-semicanonicalize before passing to libtrans.");
        //wfn_->semicanonicalize();
        Cb_ = wfn_->Cb();
    }
    process_eigenvectors();

    // Return DPD control to the user
    dpd_set_default(currentActiveDPD);

    initialized_ = true;
}

BiortIntTransform::~BiortIntTransform()
{
    if (initialized_) {
        dpd_close(myDPDNum_);
        free_int_matrix(cacheList_);
        free(cacheFiles_);
        free(zeros_);
        free(aQT_);
        free(aCorrToPitzer_);
        if(transformationType_ != TransformationType::Restricted){
            free(bQT_);
            free(bCorrToPitzer_);
        }
    }
    if(tpdm_buffer_)
        delete [] tpdm_buffer_;
}

}} // End namespaces

