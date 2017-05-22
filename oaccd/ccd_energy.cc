/*
 * @BEGIN LICENSE
 *
 * oaccd by Rolf H. Myhre May 2017, a plugin to:
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

/*Calculate the CCD Omega2 amlplitudes and return the energy.
 * Do not assume diagonal Fock matrix or orthogonal orbitals.
 * Basically use the standard CCSD algorithm for T1-transformed integals
 * PS. Use the naming convention from Helgaker et al. for the different terms*/

#include "oaccd.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"

using namespace std;

namespace psi{ namespace oaccd {


double Oaccd::ccd_energy_rhf(){

    outfile->Printf("Hej hej hallo\n");

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    psio_->open(PSIF_CC_TAMPS, PSIO_OPEN_OLD);
    psio_->open(PSIF_CC_HBAR, PSIO_OPEN_NEW);

    dpdbuf4 T2; //old T2 amplitudes
    dpdbuf4 Omega2; //Omega2 amplitudes

    ccd_a2_rhf();

    ccd_b2_rhf();

    ccd_c2_rhf();

    ccd_d2_rhf();

    ccd_e2_rhf();

    psio_->close(PSIF_CC_HBAR,true); //full of garbage
    psio_->close(PSIF_CC_TAMPS,true);
    psio_->close(PSIF_LIBTRANS_DPD,true);

//    global_dpd_->buf4_print(&T2old, "outfile", 1);
    return 0.0;

}

void Oaccd::ccd_a2_rhf(){

    //A2 term: Omega^ab_ij += (ai|bj) + Sum_cd t^cd_ij*(ac|bd), the most expensive 

    dpdbuf4 I;  //Integrals
    dpdbuf4 T2; //old T2 amplitudes
    dpdbuf4 Omega2; //Omega2 vector

    timer_on("A2");

    //Open a buffer for the aibj integrals
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "(VO|VO) (i,j,a,b)");

    //Copy integrals into amplitude buffer, first term
    global_dpd_->buf4_copy(&I, PSIF_CC_TAMPS, "Omega2");

    global_dpd_->buf4_close(&I);

    //Open a buffer for the old amplitudes
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2");

    //Open a buffer for the new amplitudes
    global_dpd_->buf4_init(&Omega2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Omega2");

    //Open a buffer for the abcd integrals
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "(VV|VV) (a,c,b,d)");

    //Contract old amplitudes with integrals
    global_dpd_->contract444(&T2,&I,&Omega2,0,0,1.0,1.0);

    global_dpd_->buf4_sort(&Omega2, PSIF_CC_TAMPS , prqs, ID("[O,V]"), ID("[O,V]"), "Omega2 iajb");

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&Omega2);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_init(&Omega2, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Omega2 iajb");
    outfile->Printf("\n\nAfter A lalala\n");
    global_dpd_->buf4_print(&Omega2, "outfile", 1);
    global_dpd_->buf4_close(&Omega2);

    timer_off("A2");
    
}

void Oaccd::ccd_b2_rhf(){

    //B2 term: Omega^ab_ij += Sum_kl t^ab_kl((ki|lj) + Sum_cd t^cd_ij*(kc|ld))

    dpdbuf4 I;      //Integrals, mostly
    dpdbuf4 W;      //Intermediates, mostly
    dpdbuf4 Omega2; //Omega2
    dpdbuf4 T2;     //T2 amplitudes

    timer_on("B2");

    //Open a buffer for the old amplitudes
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2");

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "(OV|OV) (i,j,a,b)");

    //Open an intermediate buffer
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "W OOOO");

    //Contract old amplitudes with integrals
    //W is stored as (i,j,k,l) where k and l will be contracted 
    //Sorry forinconsistent notation
    global_dpd_->contract444(&T2,&I,&W,0,0,1.0,0.0);
    global_dpd_->buf4_close(&I);

    //Read in (ij|kl) integrals and add
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "(OO|OO) (i,k,j,l)");

    global_dpd_->buf4_axpy(&I, &W, 1.0);
    global_dpd_->buf4_close(&I);

    //Contract again t^ab_ij += Sum_kl t^ab_kl*W_kilj
    global_dpd_->buf4_init(&Omega2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Omega2");

    global_dpd_->contract444(&W,&T2,&Omega2,0,1,1.0,1.0);

    global_dpd_->buf4_sort(&Omega2, PSIF_CC_TAMPS , prqs, ID("[O,V]"), ID("[O,V]"), "Omega2 iajb");

    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Omega2);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_init(&Omega2, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Omega2 iajb");
    outfile->Printf("\n\nAfter B lalala\n");
    global_dpd_->buf4_print(&Omega2, "outfile", 1);
    global_dpd_->buf4_close(&Omega2);

    timer_off("B2");
    
}

void Oaccd::ccd_c2_rhf(){

    //C2 term: Omega^ab_ij += -0.5*P^ab_ij Sum_ck t^bc_kj((ki|ac) - 0.5*Sum_dl t^ad_lj*(kd|lc))

    dpdbuf4 I;      //Integrals, mostly
    dpdbuf4 W;      //Intermediates, mostly
    dpdbuf4 Omega2; //Omega2
    dpdbuf4 T2;     //T2 amplitudes

    timer_on("C2");

    //Resort the T^ab_ij from (i,j,a,b) to (j,a,i,b)
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2");
    global_dpd_->buf4_sort(&T2, PSIF_CC_TAMPS , qrps, ID("[O,V]"), ID("[O,V]"), "T2");
    global_dpd_->buf4_close(&T2);



    //Sort (OV|OV) integrals from (i,j,a,b) to (i,b,j,a) order
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "(OV|OV) (i,j,a,b)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, psqr, ID("[O,V]"), ID("[O,V]"), "(OV|OV) (i,b,j,a)");
    global_dpd_->buf4_close(&I);



    //Open an intermediate buffer and amplitudes
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "(OV|OV) (i,b,j,a)");


    //Contract amplitudes with integrals
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W OVOV");
    global_dpd_->contract444(&T2, &I, &W, 0, 1, -0.5, 0.0);
    global_dpd_->buf4_close(&I);


    //Add (ij|ab) integrals, 
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "(VV|OO) (j,a,i,b)");
    global_dpd_->buf4_axpy(&I, &W, 1.0);
    global_dpd_->buf4_close(&I);


    //Contract with amplitudes
    global_dpd_->buf4_init(&I, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W OVOV2");
    global_dpd_->contract444(&W, &T2, &I, 0, 0, -0.5, 0.0);
    global_dpd_->buf4_close(&W);


    //Switch i and j and add
    global_dpd_->buf4_sort(&I, PSIF_CC_HBAR, rqps, ID("[O,V]"), ID("[O,V]"), "W OVOV");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W OVOV");
    global_dpd_->buf4_axpy(&W, &I, 2.0);
    global_dpd_->buf4_close(&W);


    //Finally, do the ai, bj permutation
    global_dpd_->buf4_sort(&I, PSIF_CC_HBAR, rspq, ID("[O,V]"), ID("[O,V]"), "W OVOV");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W OVOV");
    global_dpd_->buf4_axpy(&W, &I, 1.0);
    global_dpd_->buf4_close(&W);


    //Add to Omega
    global_dpd_->buf4_init(&Omega2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Omega2");
    global_dpd_->buf4_sort(&Omega2, PSIF_CC_TAMPS, prqs, ID("[O,V]"), ID("[O,V]"), "Omega2");
    global_dpd_->buf4_close(&Omega2);
    global_dpd_->buf4_init(&Omega2, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Omega2");
    global_dpd_->buf4_axpy(&I, &Omega2, 1.0);


    outfile->Printf("\n\nAfter C lalala\n");
    global_dpd_->buf4_print(&Omega2, "outfile", 1);

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&Omega2);


    timer_off("C2");
    

}

void Oaccd::ccd_d2_rhf(){

    //D2 term: t^ab_ij += 0.5*P^ab_ij Sum_ck u^bc_jk(L_aikc + 0.5*Sum_dl u^ad_li*L_ldkc)

    dpdbuf4 I;      //Integrals, mostly
    dpdbuf4 W;      //Intermediates, mostly
    dpdbuf4 Omega2; //Omega2
    dpdbuf4 T2;     //T2 amplitudes

    timer_on("D2");

    //Construct L_iajb = 2(ia|jb) - (ib|ja)
    //(OV|OV) sorted to i,b,j,a order in c term 
    global_dpd_->buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "(OV|OV) (i,b,j,a)");
    global_dpd_->buf4_sort(&W, PSIF_LIBTRANS_DPD, psrq, ID("[O,V]"), ID("[O,V]"), "L_pqrs");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "L_pqrs");
    global_dpd_->buf4_scm(&I, 2.0);
    global_dpd_->buf4_axpy(&W, &I, -1.0);
    global_dpd_->buf4_close(&W);

    //Construct u^ab_ij = 2t^ab_ij - t^ba_ij
    //T2 amplitudes in j,a,i,b order after C2 term
    global_dpd_->buf4_init(&W, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2");
    global_dpd_->buf4_sort(&W, PSIF_CC_TAMPS, rqps, ID("[O,V]"), ID("[O,V]"), "u2");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "u2");
    global_dpd_->buf4_scm(&T2, 2.0);
    global_dpd_->buf4_axpy(&W, &T2, -1.0);
    global_dpd_->buf4_close(&W);

    //Contract u and L
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W OVOV");
    global_dpd_->contract444(&T2, &I, &W, 0, 1, 0.5, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&T2);


    //Construct L_aijb = 2(ai|jb) - (ab|ji)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "(VV|OO) (j,a,i,b)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rqps, ID("[O,V]"), ID("[O,V]"), "(VV|OO) (i,a,j,b)");
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "(VO|OV) (i,a,j,b)");
    global_dpd_->buf4_scm(&I, 2.0);
    global_dpd_->buf4_init(&T2, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "(VV|OO) (i,a,j,b)");
    global_dpd_->buf4_axpy(&T2, &I, -1.0);
    global_dpd_->buf4_close(&T2);

    //Add to intermediate
    global_dpd_->buf4_axpy(&I, &W, 1.0);
    global_dpd_->buf4_close(&I);

    //Contract intermediate and u
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "u2");
    global_dpd_->buf4_init(&I, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W OVOV2");
    global_dpd_->contract444(&W, &T2, &I, 0, 0, 0.5, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    //Permute abij and add
    global_dpd_->buf4_sort(&I, PSIF_CC_HBAR, rspq, ID("[O,V]"), ID("[O,V]"), "W OVOV");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W OVOV");
    global_dpd_->buf4_axpy(&I, &W, 1.0);
    global_dpd_->buf4_close(&I);

    //Add to Omega2
    global_dpd_->buf4_init(&Omega2, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Omega2");
    global_dpd_->buf4_axpy(&W, &Omega2, 1.0);

    outfile->Printf("\n\nAfter D lalala\n");
    global_dpd_->buf4_print(&Omega2, "outfile", 1);

    global_dpd_->buf4_close(&Omega2);
    global_dpd_->buf4_close(&W);

    timer_off("D2");

}

void Oaccd::ccd_e2_rhf(){


    //E2 term Omega2 += P^ab_ij sum_c t^ac_ij*(F_bc - sum_dkl u^bd_kl*(ld|kc)) 
    //                - P^ab_ij sum_k t^ab_ik*(F_kj + sum_cdl u^cd_lj*(kd|lc))

    dpdbuf4 I;      //Integrals, mostly
    dpdbuf4 W;      //Intermediates, mostly
    dpdbuf4 Omega2; //Omega2
    dpdbuf4 T2;     //T2 amplitudes
    dpdfile2 F,D;   //Two index stuff

    timer_on("E2");
    
    //Do u^cd_lj*(kd|lc) first
    //Open u^ab_ij (i,a,j,b) and (ia|jv) (i,b,j,a) and contract
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "u2");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "(OV|OV) (i,b,j,a)");


    global_dpd_->file2_init(&D, PSIF_CC_HBAR, 0, ID('O'), ID('O'), "W OO");
    global_dpd_->contract442(&I, &T2, &D, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&T2);


    //Get Occupie-Occupied Fock block and add
    global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F OO");
    global_dpd_->file2_axpy(&F, &D, 1.0, false);
    global_dpd_->file2_close(&F);


    //Get T2 amplitudes and sort them to i,a,j,b order
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2");
    global_dpd_->buf4_sort(&T2, PSIF_CC_TAMPS, rqps, ID("[O,V]"), ID("[O,V]"), "T2");


    //Contract T2 and intermediate
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W OVOV");
    global_dpd_->contract424(&T2, &D, &W, 3, 0, false, -1.0, 0.0);

    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&D);


    //Get amplitudes and integrals and sort integrals
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "u2");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "(OV|OV) (i,b,j,a)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, psrq, ID("[O,V]"), ID("[O,V]"), "(OV|OV) (i,b,j,a)");


    //Contract integrals to intermediate
    global_dpd_->file2_init(&D, PSIF_CC_HBAR, 0, ID('V'), ID('V'), "W VV");
    global_dpd_->contract442(&T2, &I, &D, 1, 3, -1.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&T2);


    //Get virtual-virtual Fock block and add
    global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F VV");
    global_dpd_->file2_axpy(&F, &D, 1.0, false);
    global_dpd_->file2_close(&F);


    //Get T2 amplitudes and contract
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W OVOV");
    global_dpd_->contract424(&T2, &D, &W, 3, 1, false, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&D);


    //Add permutations
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, rspq, ID("[O,V]"), ID("[O,V]"), "W OVOV2");
    global_dpd_->buf4_init(&I, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W OVOV2");
    global_dpd_->buf4_axpy(&I, &W, 1.0);
    global_dpd_->buf4_close(&I);


    //Add to Omega2
    global_dpd_->buf4_init(&Omega2, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Omega2");
    global_dpd_->buf4_axpy(&W, &Omega2, 1.0);

    outfile->Printf("\n\nAfter E lalala\n");
    global_dpd_->buf4_print(&Omega2, "outfile", 1);

    global_dpd_->buf4_close(&Omega2);
    global_dpd_->buf4_close(&W);


    timer_off("E2");

}

}}//Namespaces
