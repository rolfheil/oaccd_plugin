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

#ifndef oaccd_h
#define oaccd_h

#include "psi4/psi4-dec.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"

using namespace std;

namespace psi{ namespace oaccd {

class Oaccd : public Wavefunction
{
public:
    Oaccd(SharedWavefunction ref_wfn, Options& options);
    virtual ~Oaccd();

    double compute_energy();

private:
    Dimension virtpi_;
    void common_init();

    void int_trans_rhf();

    //Various option variables
    string reference;

    //Library stuff
    IntegralTransform *ints;

    //Tensors
    SharedMatrix MOoeIntsA;
    SharedMatrix FockA;

};


}} // End namespaces

#endif //oaccd.h
