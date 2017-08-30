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

    common_initialize();

    if(init) initialize();
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
        if(transformationType_ != Restricted){
            free(bQT_);
            free(bCorrToPitzer_);
        }
    }
    if(tpdm_buffer_)
        delete [] tpdm_buffer_;
}

}} // End namespaces

