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

#include "biortwfn.h"

using namespace std;

namespace psi{ namespace biortwfn {

//Biortwfn::Biortwfn(SharedWavefunction ref_wfn, Options& options)
//        : Wavefunction(options)
//{   
//    // Shallow copy ref_wfn data into this wavefunction                                 
//    shallow_copy(ref_wfn);
//}
  

SharedMatrix Biortwfn::lCa() const {
    if (!lCa_) {
        if (!reference_wavefunction_)
            throw PSIEXCEPTION("Biortwfn::lCa: Unable to obtain MO coefficients.");
        else
            return reference_wavefunction_->Ca();
    }

    return lCa_;
}

SharedMatrix Biortwfn::rCa() const {
    if (!rCa_) {
        if (!reference_wavefunction_)
            throw PSIEXCEPTION("Biortwfn::rCa: Unable to obtain MO coefficients.");
        else
            return reference_wavefunction_->Ca();
    }

    return rCa_;
}

SharedMatrix Biortwfn::lCb() const {
    if (!lCb_) {
        if (!reference_wavefunction_)
            throw PSIEXCEPTION("Biortwfn::lCb: Unable to obtain MO coefficients.");
        else
            return reference_wavefunction_->Ca();
    }

    return lCb_;
}

SharedMatrix Biortwfn::rCb() const {
    if (!rCb_) {
        if (!reference_wavefunction_)
            throw PSIEXCEPTION("Biortwfn::rCb: Unable to obtain MO coefficients.");
        else
            return reference_wavefunction_->Ca();
    }

    return rCb_;
}

}} // End namespaces

