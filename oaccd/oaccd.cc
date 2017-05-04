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
#include "psi4/libdpd/dpd.h"

using namespace std;

namespace psi{ namespace oaccd {

//Will eventually become NOCCD and OACCD routine
//Currently working on the CCD part
//This is similar to the CCSD and OCC stuff in Psi4

Oaccd::Oaccd(SharedWavefunction ref_wfn, Options& options)
    : Wavefunction(options)
{
    // Shallow copy ref_wfn data into this wavefunction
    shallow_copy(ref_wfn);
    common_init();
}

Oaccd::~Oaccd()
{
}

void Oaccd::common_init()
{
    // nsopi_, frzcpi_, etc are Dimension objects for symmetry orbitals
    // These are copied from ref_wfn when we call for shallow_copy
    virtpi_ = nsopi_ - frzcpi_ - frzvpi_ - doccpi_;
    adoccpi_ = doccpi_ - frzcpi_;
    avirtpi_ = virtpi_ - frzvpi_;

    outfile->Printf("The wavefunction has the following dimensions:\n");
    nsopi_.print();
    frzcpi_.print();
    doccpi_.print();
    virtpi_.print();
    frzvpi_.print();
    nmopi_.print();
    outfile->Printf("Irreps: %3i \n",nirrep_);

    reference=options_.get_str("REFERENCE");

    if(reference == "RHF") {//Only RHF for the time being

        //Allocate memory
        MOoeIntsA = std::shared_ptr<Matrix>(
        new Matrix("MO-basis alpha one-electron ints", nirrep_, nmopi_, nmopi_));
        FockA = std::shared_ptr<Matrix>(
        new Matrix("MO-basis alpha Fock matrix", nirrep_, nmopi_, nmopi_));

        FDiaOccA = std::shared_ptr<Vector>(        
        new Vector("Fock matrix occupied diagonal", adoccpi_));
        FDiaVirA = std::shared_ptr<Vector>(        
        new Vector("Fock matrix virtual diagonal", avirtpi_));

//        FockA = Fa();
        Fa_->print();
        Ca_->print();
    }
    else{
        throw PSIEXCEPTION("OACCD only implemented for RHF");
    }


}

double Oaccd::compute_energy()
{
    //Start by getting the required integrals

    //Allocate integrals,must be done after constructor
    std::vector<std::shared_ptr<MOSpace>> spaces = {MOSpace::occ, MOSpace::vir};
    ints = new IntegralTransform(shared_from_this(), spaces,
               IntegralTransform::Restricted,
               IntegralTransform::DPDOnly,
               IntegralTransform::QTOrder,
               IntegralTransform::None,
               false);

    ints->set_dpd_id(0);    
    ints->set_keep_dpd_so_ints(true);
    ints->initialize();
    dpd_set_default(ints->get_dpd_id());

    int_trans_rhf();
                   
//    Fa_->print();
//    FockA = Fa_;
//    FockA->transform(Ca_);
//    FockA->print();

    //Generete a start guess. MP2? Might be bad for difficult cases
    

    //Start orbital iteration loop here 
    

    //Start DIIS CC iterations. Implement amplitude equations and figure out DIIS library

    
    //Orbital gradients
    

    //Check convergence, if not, update for next iteration
    

    //End orbital loop    

    return 0.0;
}

extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "OACCD"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C"
SharedWavefunction oaccd(SharedWavefunction ref_wfn, Options& options)
{

    // Note that if this function was integrated into Psi4 we would not be using P::e.wavefunction
    // Instead everything would be explicitly passed
    SharedWavefunction wfn(new Oaccd(ref_wfn, options));
    wfn->compute_energy();

    return wfn;
}

}} // End namespaces

