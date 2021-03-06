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

#ifndef OACCD_H
#define OACCD_H

#include "defines.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/local.h"
#include "psi4/libdiis/diismanager.h"
#include "biort_libtrans/biort_inttransform.h"
#include "biort_wfn/biortwfn.h"

using namespace std;

namespace psi{ namespace oaccd {

class Oaccd : public biortwfn::Biortwfn
{
public:
    Oaccd(SharedWavefunction ref_wfn, Options& options);
    ~Oaccd();

    double compute_energy();

protected:

    void common_init();

    void int_trans_rhf();
    void f_denominator();
    double mp2_energy_rhf();

    double ccd_energy_rhf();
    void ccd_a2_rhf();
    void ccd_b2_rhf();
    void ccd_c2_rhf();
    void ccd_d2_rhf();
    void ccd_e2_rhf();
    void ccd_L_rhf();
    double ccd_update_rhf();

    //Energies etc.
    double oaccd_energy;
    double tF_energy; //Transformed reference energy
    double nuc_energy; //Nuclear repulsion energy
    double mp2_energy;
    double ccd_energy;
    double omega_norm;
    double temp_norm;

    //Various option variables
    double r_convergence;
    double e_convergence;
    int cc_maxiter;
    int cc_maxdiis;
    int cc_mindiis;
    string reference;
    int print_;

    //Library stuff
    BiortIntTransform *ints;
    DIISManager *t2DiisManager;

    //Tensors
    SharedMatrix MOoeIntsA;
    SharedMatrix FockA;

    //Fock diagonals
    SharedVector FDiaOccA;
    SharedVector FDiaVirA;

    //Dimension stuff
    Dimension virtpi_;
    Dimension adoccpi_;
    Dimension avirtpi_;
};


}} // End namespaces

#endif //oaccd.h
