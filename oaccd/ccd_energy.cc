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

    //Sort amplitudes to iajb order
    global_dpd_->buf4_init(&Omega2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Omega ijab");
    global_dpd_->buf4_sort(&Omega2, PSIF_CC_TAMPS , prqs, ID("[O,V]"), ID("[O,V]"), "Omega iajb");
    global_dpd_->buf4_close(&Omega2);


    //Sort the old amplitudes
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T ijab");
    global_dpd_->buf4_sort(&T2, PSIF_CC_TAMPS , prqs, ID("[O,V]"), ID("[O,V]"), "T iajb");
    global_dpd_->buf4_close(&T2);

    ccd_c2_rhf();

    ccd_d2_rhf();

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
                  ID("[O,O]"), ID("[V,V]"), 0, "g_aibj <OO|VV>");

    //Copy integrals into amplitude buffer, first term
    global_dpd_->buf4_copy(&I, PSIF_CC_TAMPS, "Omega ijab");

    global_dpd_->buf4_close(&I);

    //Open a buffer for the old amplitudes
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T ijab");

    //Open a buffer for the new amplitudes
    global_dpd_->buf4_init(&Omega2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Omega ijab");

    //Open a buffer for the abcd integrals
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "g_acbd <VV|VV>");

    //Contract old amplitudes with integrals
    global_dpd_->contract444(&T2,&I,&Omega2,0,0,1.0,1.0);

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&Omega2);
    global_dpd_->buf4_close(&T2);

    timer_off("A2");
    
}

void Oaccd::ccd_b2_rhf(){

    //B2 term: Omega^ab_ij += Sum_kl t^ab_kl((ki|lj) + Sum_cd t^cd_ij*(kc|ld))

    dpdbuf4 I;      //Integrals
    dpdbuf4 W;      //Intermediate
    dpdbuf4 Omega2; //Omega2 amplitudes
    dpdbuf4 T2;     //old T2 amplitudes

    timer_on("B2");

    //Open a buffer for the old amplitudes
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T ijab");

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "g_iajb <OO|VV>");

    //Open an intermediate buffer
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "W ijkl");

    //Contract old amplitudes with integrals
    global_dpd_->contract444(&T2,&I,&W,0,0,1.0,0.0);

    global_dpd_->buf4_close(&I);

    //Read in g_ikjl integrals and add
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "g_ikjl <OO|OO>");

    global_dpd_->buf4_axpy(&I, &W, 1.0);
    global_dpd_->buf4_close(&I);

    //Contract again t^ab_ij += Sum_kl t^ab_kl*W_kilj
    global_dpd_->buf4_init(&Omega2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Omega ijab");

    global_dpd_->contract444(&W,&T2,&Omega2,1,1,1.0,1.0);

    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Omega2);
    global_dpd_->buf4_close(&T2);

    timer_off("B2");
    
}

void Oaccd::ccd_c2_rhf(){

    //C2a term: Omega^ab_ij += -0.5*P^ab_ij Sum_ck t^bc_kj((ki|ac) - 0.5*Sum_dl t^ad_lj*(kd|lc))

    dpdbuf4 I;      //Integrals
    dpdbuf4 W;      //Intermediate
    dpdbuf4 Omega2; //Omega2 amplitudes
    dpdbuf4 T2;     //T2 amplitudes

    timer_on("C2");
    timer_on("C2a");

    //Open an intermediate buffer and old amplitudes
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T iajb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W iajb");

    //Get the (ia|jb) integrals, use this file address because it's cheaper to 
    //sort one set to another than transforming an extra set
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "g_iajb <OO|VV>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "g_iajb (OV|OV)");
    global_dpd_->buf4_close(&I);

    //Contract old amplitudes with integrals
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "g_iajb (OV|OV)");
    global_dpd_->contract444(&T2, &I, &W, 0, 0, -0.5, 0.0);
    global_dpd_->buf4_close(&I);

    //Add (ki|ac) integrals
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "g_ijab (OV|OV)");
    global_dpd_->buf4_axpy(&I, &W, 1.0);
    global_dpd_->buf4_close(&I);

    
    //Contract with old amplitudes
    global_dpd_->buf4_init(&I, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W C2a");
    global_dpd_->contract444(&W, &T2, &I, 0, 0, -0.5, 0.0);
    global_dpd_->buf4_close(&W);

    timer_off("C2a");
    
    //C2b term: Omega^ab_ij += -P^ab_ij Sum_ck t^bc_ki((kj|ac) - 0.5*Sum_dl t^ad_lj*(kd|lc))

    timer_on("C2b");

    //Open an intermediate buffer and old amplitudes
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W iajb");

    //Get the iajb integrals. They were sorted nthe C term
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "g_iajb (OV|OV)");
    global_dpd_->contract444(&T2, &I, &W, 0, 0, -0.5, 0.0);
    global_dpd_->buf4_close(&I);

    //Add (kj|ac) integrals
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "g_ijab (OV|OV)");
    global_dpd_->buf4_axpy(&I, &W, 1.0);
    global_dpd_->buf4_close(&I);

    //Contract with amplitudes
    global_dpd_->buf4_init(&I, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W C2b");
    global_dpd_->contract444(&W, &T2, &I,0,0,-0.5,0.0);
    global_dpd_->buf4_close(&W);

    //Permute ij
    global_dpd_->buf4_sort(&I, PSIF_CC_HBAR , rqps, ID("[O,V]"), ID("[O,V]"), "W iajb");
    global_dpd_->buf4_close(&I);

    timer_off("C2b");

    //Add C2a and C2b
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W C2a");
    global_dpd_->buf4_init(&I, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W iajb");
    global_dpd_->buf4_axpy(&I, &W, 1.0);
    global_dpd_->buf4_close(&I);

    //Get permutation
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR , rspq, ID("[O,V]"), ID("[O,V]"), "W iajb");
    global_dpd_->buf4_init(&I, PSIF_CC_HBAR, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W iajb");
    global_dpd_->buf4_axpy(&I, &W, 1.0);

    global_dpd_->buf4_init(&Omega2, PSIF_CC_TAMPS, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Omega iajb");
    global_dpd_->buf4_axpy(&W, &Omega2, -1.0);

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&Omega2);
    global_dpd_->buf4_close(&T2);

    timer_off("C2");

}

void Oaccd::ccd_d2_rhf(){

    //D2 term: t^ab_ij += 0.5*P^ab_ij Sum_ck u^bc_jk(L_aikc + 0.5*Sum_dl u^ad_li*L_ldkc)

    dpdbuf4 I;      //Integrals
    dpdbuf4 W;      //Intermediate
    dpdbuf4 T2;     //T2 amplitudes
    dpdbuf4 T2old;  //old T2 amplitudes

}

}}//Namespaces
