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

using namespace std;

namespace psi{ namespace oaccd {

void Oaccd::int_trans_rhf(){

//Routine to update integrals, Fock matrix etc. to the current basis set
//Might be necessary to inherit a new IntegralTransform class to account for the 
//biorthogonal transformation

    outfile->Printf("Hej hej hallaa\n");
    
    //Update one-electron orbitals etc. 
    //Coefficients from Ca_ (and Cb_) in wafunctionobject
    
    ints->set_print(5);
    ints->update_orbitals();
    ints->set_keep_dpd_so_ints(1);

    //Transform to (OO|OO), must figure out how to do it when left and right is different
    timer_on("Trans (OO|OO)");
    ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ, 
                        IntegralTransform::MakeAndKeep);
    timer_off("Trans (OO|OO)");

    //Transform to (OO|VV)
    timer_on("Trans (OO|VV)");
    ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir, 
                        IntegralTransform::ReadAndNuke);
    timer_off("Trans (OO|VV)");

    //Transform to (OV|OV)
    timer_on("Trans (OV|OV)");
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir, 
                        IntegralTransform::MakeAndNuke);
    timer_off("Trans (OV|OV)");

    //Transform to (VO|VO) not the same as above in NOCC
    timer_on("Trans (VO|VO)");
    ints->transform_tei(MOSpace::vir, MOSpace::occ, MOSpace::vir, MOSpace::occ, 
                        IntegralTransform::MakeAndNuke);
    timer_off("Trans (VO|VO)");

    //Transform to (VV|VV)
    timer_on("Trans (VV|VV)");
    ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::vir, MOSpace::vir, 
                        IntegralTransform::MakeAndNuke);
    timer_off("Trans (VV|VV)");

}

}} // End namespaces
