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

#include "oaccd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libpsio/psio.hpp"

using namespace std;

namespace psi{ namespace oaccd {

void Oaccd::int_trans_rhf(){

//Routine to update integrals, Fock matrix etc. to the current basis set
//Might be necessary to inherit a new IntegralTransform class to account for the 
//biorthogonal transformation

    dpdbuf4 K, G;

    //Update one-electron orbitals etc. 
    //Coefficients from Ca_ (and Cb_) in wafunctionobject
    
    ints->set_print(0);
    ints->update_orbitals();
    ints->set_keep_dpd_so_ints(1);

    //Transform to (OO|OO), must write a new transform class for NO transformations
    timer_on("Trans (OO|OO)");
    ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ, 
                        IntegralTransform::MakeAndNuke);
    timer_off("Trans (OO|OO)");

    //Transform to (OV|OV)
    timer_on("Trans (OV|OV)");
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir, 
                        IntegralTransform::MakeAndNuke);
    timer_off("Trans (OV|OV)");

    //Transform to (OO|VV)
    timer_on("Trans (VV|OO)");
    ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::occ, 
                        IntegralTransform::MakeAndKeep);
    timer_off("Trans (VV|OO)");

    //Transform to (VV|VV)
    timer_on("Trans (VV|VV)");
    ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::vir, MOSpace::vir, 
                        IntegralTransform::ReadAndNuke);
    timer_off("Trans (VV|VV)");

    //DPD needs both contracted indices in either row or column, so sort to 
    //physics notation

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // (OO|OO) -> <OO|OO>
    timer_on("Sort (OO|OO) -> <OO|OO>");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                 ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[O,O]"), "g_ikjl <OO|OO>");
    global_dpd_->buf4_close(&K);
    timer_off("Sort (OO|OO) -> <OO|OO>");

    // (VV|VV) -> <VV|VV>
    timer_on("Sort (VV|VV) -> <VV|VV>");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                 ID("[V>=V]+"), ID("[V>=V]+"), 0, "MO Ints (VV|VV)");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[V,V]"), ID("[V,V]"), "g_acbd <VV|VV>");
    global_dpd_->buf4_close(&K);
    timer_off("Sort (VV|VV) -> <VV|VV>");

    // (VV|OO) -> (OO|VV)
    timer_on("Sort (VV|OO) -> (OV|OV)");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                 ID("[V>=V]+"), ID("[O>=O]+"), 0, "MO Ints (VV|OO)");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD, sprq, ID("[O,V]"), ID("[O,V]"), "g_ijab (OV|OV)");
    global_dpd_->buf4_close(&K);
    timer_off("Sort (VV|OO) -> (OV|OV)");

    // (OV|OV) -> <OO|VV>
    timer_on("Sort (OV|OV) -> <OO|VV>");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                 ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[V,V]"), "g_iajb <OO|VV>");
    global_dpd_->buf4_close(&K);
    timer_off("Sort (OV|OV) -> <OO|VV>");

    //Construct L_iajb = 2g_iajb - g_ibja
    timer_on("Construct L_iajb");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                 ID("[O,O]"), ID("[V,V]"), 0, "g_iajb <OO|VV>");

    //Copy to g_aibj <OO|VV>, they are the same in orthogonal basis
    global_dpd_->buf4_copy(&K, PSIF_LIBTRANS_DPD, "g_aibj <OO|VV>");

    global_dpd_->buf4_copy(&K, PSIF_LIBTRANS_DPD, "L_iajb <OO|VV>");
    global_dpd_->buf4_close(&K);

    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                 ID("[O,O]"), ID("[V,V]"), 0, "L_iajb <OO|VV>");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "g_ibja <OO|VV>");
    global_dpd_->buf4_scm(&K, 2.0);

    global_dpd_->buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                 ID("[O,O]"), ID("[V,V]"), 0, "g_ibja <OO|VV>");
    global_dpd_->buf4_axpy(&G, &K, -1.0);

    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&K);
    timer_off("Construct L_iajb");

    timer_on("F denominator");
    f_denominator();
    timer_off("F denominator");
    
    psio_->close(PSIF_LIBTRANS_DPD,true);
}

void Oaccd::f_denominator(){

    //Make a four index buffer
    dpdbuf4 D;

    //Transform the Fock matrix to MO basis
    FockA = Fa_;
    FockA->transform(Ca_);
    FockA->print();

    for(int h = 0; h < nirrep_; ++h){
        for(int i = frzcpi_[h]; i < doccpi_[h]; i++){
            FDiaOccA->set(h,i - frzcpi_[h],FockA->get(h,i,i));
        }
        for(int a = doccpi_[h]; a < doccpi_[h] + avirtpi_[h]; a++){
            FDiaVirA->set(h,a - doccpi_[h],FockA->get(h,a,a));
        }
    }
//    outfile->Printf("\n lalala %f \n", FDiaOccA[4]);
    FDiaOccA->print();
    FDiaVirA->print();

    //Build the denominators
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&D, h);
        for(int row = 0;row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(FDiaOccA->get(i) + FDiaOccA->get(j) - 
                                             FDiaVirA->get(a) - FDiaVirA->get(b));
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&D, h);
        global_dpd_->buf4_mat_irrep_close(&D, h);
    }
    global_dpd_->buf4_close(&D);

}
}} // End namespaces
