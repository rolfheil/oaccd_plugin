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

    //Transform to (VO|VO)
    timer_on("Trans (VO|VO)");
    ints->transform_tei(MOSpace::vir, MOSpace::occ, MOSpace::vir, MOSpace::occ, 
                        IntegralTransform::MakeAndKeep);
    timer_off("Trans (VO|VO)");

    //Transform to (VO|OV)
    timer_on("Trans (VO|OV)");
    ints->transform_tei(MOSpace::vir, MOSpace::occ, MOSpace::occ, MOSpace::vir, 
                        IntegralTransform::ReadAndNuke);
    timer_off("Trans (VO|OV)");

    //Transform to (OV|OV)
    timer_on("Trans (OV|OV)");
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir, 
                        IntegralTransform::MakeAndNuke);
    timer_off("Trans (OV|OV)");

    //Transform to (VV|OO)
    timer_on("Trans (VV|OO)");
    ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::occ, 
                        IntegralTransform::MakeAndKeep);
    timer_off("Trans (VV|OO)");

    //Transform to (VV|VV)
    timer_on("Trans (VV|VV)");
    ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::vir, MOSpace::vir, 
                        IntegralTransform::ReadAndNuke);
    timer_off("Trans (VV|VV)");

    //DPD needs both contracted indices in either row or column, so sort the integrals. 
    //In biorthogonal basis,, (VO|VO) =/= (OV|OV)

    outfile->Printf("\n tF_energy is: %f \n",tF_energy);
    tFa_ = ints->compute_biort_fock_matrix(H_,lCa_,rCa_,tF_energy);
    tF_energy = tF_energy + nuc_energy; 
    outfile->Printf("\n tF_energy is: %f \n",tF_energy);

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // (l,j,k,i) -> (i,j,k,l) 
    timer_on("Sort (OO|OO) (l,j,k,i) -> (i,j,k,l)");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                 ID("[O,O]"), ID("[O,O]"), 0, "MO Ints (OO|OO)");
    outfile->Printf("OOOO integrals");
    global_dpd_->buf4_print(&K,"tull3",1);   
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , sqrp, ID("[O,O]"), ID("[O,O]"), "(OO|OO) (i,j,k,l)");
    global_dpd_->buf4_close(&K);
    timer_off("Sort (OO|OO) (l,j,k,i) -> (i,j,k,l)");


    // (a,b,c,d) -> (a,c,b,d)
    timer_on("Sort (VV|VV) (a,b,c,d) -> (a,c,b,d)");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                 ID("[V,V]"), ID("[V,V]"), 0, "MO Ints (VV|VV)");
    outfile->Printf("VVVV integrals");
    global_dpd_->buf4_print(&K,"tull2",1);   
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[V,V]"), ID("[V,V]"), "(VV|VV) (a,c,b,d)");
    global_dpd_->buf4_close(&K);
    timer_off("Sort (VV|VV) (a,b,c,d) -> (a,c,b,d)");


    // (a,b,i,j) -> (j,a,i,b)
    timer_on("Sort (VV|OO) (a,b,i,j) -> (j,a,i,b)");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                 ID("[V,V]"), ID("[O,O]"), 0, "MO Ints (VV|OO)");
    outfile->Printf("VVOO integrals");
    global_dpd_->buf4_print(&K,"tull2",1);   
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD, sprq, ID("[O,V]"), ID("[O,V]"), "(VV|OO) (j,a,i,b)");
    global_dpd_->buf4_close(&K);
    timer_off("Sort (VV|OO) (a,b,i,j) -> (j,a,i,b)");


    //(VO|VO) not the same as (OV|OV) 
    timer_on("Sort (VO|VO) (b,j,a,i) -> (i,j,a,b)");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[V,O]"),
                 ID("[V,O]"), ID("[V,O]"), 0, "MO Ints (VO|VO)");
    
    outfile->Printf("VOVO integrals");
    global_dpd_->buf4_print(&K,"tull",1);   
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , sqrp , ID("[O,O]"), ID("[V,V]"), "(VO|VO) (i,j,a,b)");
    global_dpd_->buf4_close(&K);
    timer_off("Sort (VO|VO) (b,j,a,i) -> (i,j,a,b)");


    // (OV|OV) -> <OO|VV>
    timer_on("Sort (OV|OV) (i,a,j,b) -> (i,j,a,b)");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                 ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    outfile->Printf("OVOV integrals");
    global_dpd_->buf4_print(&K,"tull1",1);   
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[V,V]"), "(OV|OV) (i,j,a,b)");
    global_dpd_->buf4_close(&K);
    timer_off("Sort (OV|OV) (i,a,j,b) -> (i,j,a,b)");


    //We need (VO|OV) as well in L_aijb
    timer_on("Sort (VO|OV) (a,i,j,b) -> (i,a,j,b)");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                 ID("[V,O]"), ID("[O,V]"), 0, "MO Ints (VO|OV)");
    outfile->Printf("VOOV integrals");
    global_dpd_->buf4_print(&K,"tull",1);   
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , qprs, ID("[O,V]"), ID("[O,V]"), "(VO|OV) (i,a,j,b)");
    global_dpd_->buf4_close(&K);
    timer_off("Sort (VO|OV) (a,i,j,b) -> (i,a,j,b)");

    //Generate the orbital energy denominators
    timer_on("F denominator");
    f_denominator();
    timer_off("F denominator");
    

    psio_->close(PSIF_LIBTRANS_DPD,true);
}

void Oaccd::f_denominator(){

    //Buffers
    dpdbuf4 D;
    dpdfile2 F;

    //Transform the Fock matrix to MO basis
    FockA = Fa_;
    outfile->Printf("\n Original Fock matrix lalala\n ");
    Fa_->print();
    outfile->Printf("\n One electron integrals\n ");
    H_->print();
    outfile->Printf("\n Printing the Fock matrix before\n ");
    FockA->print();
    FockA->transform(lCa_, FockA, rCa_);
    outfile->Printf("\n Printing the Fock matrix after\n ");
    FockA->print();

    for(int h = 0; h < nirrep_; ++h){
        for(int i = frzcpi_[h]; i < doccpi_[h]; i++){
            FDiaOccA->set(h,i - frzcpi_[h],FockA->get(h,i,i));
        }
        for(int a = doccpi_[h]; a < doccpi_[h] + avirtpi_[h]; a++){
            FDiaVirA->set(h,a - doccpi_[h],FockA->get(h,a,a));
        }
    }

    FDiaOccA->print();
    FDiaVirA->print();

    //Build the denominators
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D (i,j,a,b)");
    
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

    //Generate the Occupied-Occupied and Virtual-Virtual blocks
    //of the Fock matrix and store as DPD. The Fock matrix is not 
    //symmetric
    
    global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F OO");
    global_dpd_->file2_mat_init(&F);

    for(int h = 0; h < nirrep_; ++h){
        for(int i = frzcpi_[h]; i < doccpi_[h]; i++){
            for(int j = frzcpi_[h]; j < doccpi_[h]; j++){
                F.matrix[h][i][j] = FockA->get(h,i,j);
            }
        }
    }

    global_dpd_->file2_mat_wrt(&F);
    global_dpd_->file2_close(&F);

    
    //Virtual block of the Fock matrix
    global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F VV");
    global_dpd_->file2_mat_init(&F);

    for(int h = 0; h < nirrep_; ++h){
        for(int a = doccpi_[h]; a < doccpi_[h] + avirtpi_[h]; a++){
            for(int b = doccpi_[h]; b < doccpi_[h] + avirtpi_[h]; b++){
                F.matrix[h][a - doccpi_[h]][b - doccpi_[h]] = FockA->get(h,a,b);
            }
        }
    }

    global_dpd_->file2_mat_wrt(&F);
    global_dpd_->file2_close(&F);

}
}} // End namespaces
