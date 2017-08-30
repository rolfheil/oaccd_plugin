/*
 * @BEGIN LICENSE
 *
 * oaccd by Psi4 Developer, a plugin to:
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

#ifndef BIORTWFN_H
#define BIORTWFN_H

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"

using namespace std;

namespace psi{ namespace biortwfn {

class Biortwfn : public Wavefunction
{
public:

   Biortwfn(SharedWavefunction ref_wfn, Options& options);

   Biortwfn(Options & options):Wavefunction(options){};
   
   //Return the alpha C matrices
    SharedMatrix lCa() const;
    SharedMatrix rCa() const;

    //Return the beta C matrices
    SharedMatrix lCb() const;
    SharedMatrix rCb() const;

protected:

    //left and right MO coeffiecent matrices alpha
    SharedMatrix lCa_;
    SharedMatrix rCa_;

    //left and right MO coeffiecent matrices beta
    SharedMatrix lCb_;
    SharedMatrix rCb_;

};


}} // End namespaces

#endif //biortwfn.h
