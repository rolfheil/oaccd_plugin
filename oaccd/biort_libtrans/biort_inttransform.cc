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

using namespace std;

namespace psi{ namespace oaccd {

BiortIntTransform::BiortIntTransform(std::shared_ptr<Wavefunction> wfn,
                                     SpaceVec spaces,
                                     TransformationType transformationType,
                                     OutputType outputType,
                                     MOOrdering moOrdering,
                                     FrozenOrbitals frozenOrbitals,
                                     bool init):
            initialized_(false),
            psio_(_default_psio_lib_),
            wfn_(wfn),
            transformationType_(transformationType),
            uniqueSpaces_(spaces),
            moOrdering_(moOrdering),
            outputType_(outputType),
            frozenOrbitals_(frozenOrbitals),
            alreadyPresorted_(false),
            dpdIntFile_(PSIF_LIBTRANS_DPD),
            aHtIntFile_(PSIF_LIBTRANS_A_HT),
            bHtIntFile_(PSIF_LIBTRANS_B_HT),
            nTriSo_(0),
            nTriMo_(0),
            nfzc_(0),
            nfzv_(0),
            spaces_(0),
            labels_(0),
            tolerance_(1.0E-16),
            moIntFileAA_(0),
            moIntFileAB_(0),
            moIntFileBB_(0),
            myDPDNum_(1),
            print_(1),
            zeros_(0),
            sosym_(0),
            mosym_(0),
            aQT_(0),
            bQT_(0),
            cacheFiles_(0),
            cacheList_(0),
            lCa_(wfn->lCa()),
            rCa_(wfn->rCa()),
            lCb_(wfn->lCb()),
            rCb_(wfn->rCb()),
            H_(wfn->H()),
            keepIwlSoInts_(false),
            keepIwlMoTpdm_(true),
            keepDpdSoInts_(false),
            keepDpdMoTpdm_(true),
            keepHtInts_(true),
            keepHtTpdm_(true),
            tpdmAlreadyPresorted_(false),
            soIntTEIFile_(PSIF_SO_TEI)
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

}} // End namespaces

