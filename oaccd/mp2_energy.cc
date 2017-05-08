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


#include "oaccd.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libpsio/psio.hpp"

using namespace std;

namespace psi{ namespace oaccd {

/*Calculate the MP2 amplitudes and energy as a start guess
 * Assume a diagonal Fock matrix etc. for the start guess.
 * T2_{aijb} = (ia|jb)*1/(F_ii + F_jj - F_aa - F_bb) 
 * E_MP2 = sum_{iajb}T2_{aibj}*(ia|jb) */

double Oaccd::mp2_energy_rhf(){

    dpdbuf4 D;  //Fock diagonal denominator and integrals
    dpdbuf4 T2; //T2 amplitudes
    
    //Open integrals and amplitude files
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    psio_->open(PSIF_CC_TAMPS, PSIO_OPEN_OLD);

    //Get the integrals and copy into T2 buffer
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "g_iajb <OO|VV>");
    global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "Tijab (old)");
    global_dpd_->buf4_close(&D);

    //Get integrals and denominators
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0,  "Tijab (old)");
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0,  "D <OO|VV>");
    //Direct product
    global_dpd_->buf4_dirprd(&D, &T2);

    //MP2 amplitudes now in buffer
    global_dpd_->buf4_close(&D);

    //Energy requires L_iajb = 2g_iajb - g_ibja
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0,  "L_iajb <OO|VV>");
    double t2_energy = global_dpd_->buf4_dot(&T2, &D);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&T2);

    psio_->close(PSIF_CC_TAMPS,true);
    psio_->close(PSIF_LIBTRANS_DPD,true);

    return t2_energy;
}
}} // End namespaces
