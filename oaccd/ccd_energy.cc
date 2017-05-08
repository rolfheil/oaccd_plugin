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

/*Calculate the CCD T2 amlplitudes and return the energy.
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

    dpdbuf4 I;  //Integrals
    dpdbuf4 T2old; //old T2 amplitudes
    dpdbuf4 T2; //T2 amplitudes
    dpdbuf4 W; //Intermediates

    //A2 term: t^ab_ij += (ai|bj) + Sum_cd t^cd_ij*(ac|bd), the most expensive 

    //Open a buffer for the aibj integrals
    timer_on("A2");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "g_aibj <OO|VV>");

    //Copy integrals into amplitude buffer, first term
    global_dpd_->buf4_copy(&I, PSIF_CC_TAMPS, "Tijab");

    global_dpd_->buf4_close(&I);

    //Open a buffer for the old amplitudes
    global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tijab (old)");

    //Open a buffer for the new amplitudes
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tijab");

    //Open a buffer for the abcd integrals
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "g_acbd <VV|VV>");

    //Contract old amplitudes with integrals
    global_dpd_->contract444(&T2old,&I,&T2,0,0,1.0,1.0);

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&T2);
    timer_off("A2");
    
    //B2 term: t^ab_ij += Sum_kl t^ab_kl((ki|lj) + Sum_cd t^cd_ij*(kc|ld))

    outfile->Printf("Hej hej hallo\n");

    timer_on("B2");

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "g_iajb <OO|VV>");

    //Open an intermediate buffer
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "W temps");

    //Contract old amplitudes with integrals
    global_dpd_->contract444(&T2old,&I,&W,0,0,1.0,0.0);

    global_dpd_->buf4_close(&I);

    //Read in g_ikjl integrals and add
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "g_ikjl <OO|OO>");
    global_dpd_->buf4_axpy(&I, &W, 1.0);
    global_dpd_->buf4_close(&I);

    //Contract again t^ab_ij += Sum_kl t^ab_kl*W_kilj
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tijab");
    global_dpd_->contract444(&W,&T2old,&T2,1,1,1.0,1.0);

    timer_off("B2");

    outfile->Printf("Hej hej hallo\n");

    psio_->close(PSIF_CC_HBAR,false);
    psio_->close(PSIF_CC_TAMPS,true);
    psio_->close(PSIF_LIBTRANS_DPD,true);


    return 0.0;

}

}}//Namespaces
