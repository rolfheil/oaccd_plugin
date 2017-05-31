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

#define PI 3.14159265/4

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
    cc_maxdiis = options_.get_int("CC_DIIS_MAX_VECS");
    cc_mindiis = options_.get_int("CC_DIIS_MIN_VECS");

    r_convergence = options_.get_double("R_CONVERGENCE");
    e_convergence = options_.get_double("E_CONVERGENCE");

    cc_maxiter = options_.get_int("CC_MAXITER");

    outfile->Printf("\nmindiis %3i \n",cc_mindiis);
    outfile->Printf("maxdiis %3i \n",cc_maxdiis);
    outfile->Printf("maxiter %3i \n",cc_maxiter);

    outfile->Printf("\ne_convergence: %16.10f \n", e_convergence);
    outfile->Printf("r_convergence: %16.10f \n", r_convergence);

    if(reference == "RHF") {//Only RHF for the time being

        //Allocate memory
        MOoeIntsA = std::shared_ptr<Matrix>(
        new Matrix("MO-basis alpha one-electron ints", nirrep_, nmopi_, nmopi_));
        FockA = std::shared_ptr<Matrix>(
        new Matrix("MO-basis alpha Fock matrix", nirrep_, nmopi_, nmopi_));

/*        Rotation = std::shared_ptr<Matrix>(
        new Matrix("Hacky rotation matrix", nirrep_, nmopi_, nmopi_));

        XCa = std::shared_ptr<Matrix>(
        new Matrix("Hacky coefficient matrix", nirrep_, nmopi_, nmopi_));

        for(int i =0; i < nmopi_[0]; ++i){
            for(int j =0; j < nmopi_[0]; ++j){

                if(i == j){ 
                    if(i == doccpi_[0]-2 || i == doccpi_[0]-1 || i == doccpi_[0] || i == doccpi_[0]+1){
                        Rotation->set(0,i,j,cos(PI));
                    } else{
                        Rotation->set(0,i,j,1.0);
                    }
                } else if(j == i-1){
                    if(i == doccpi_[0]-1 || i == doccpi_[0]+1){
                        Rotation->set(0,i,j,sin(PI));
                    } else{
                        Rotation->set(0,i,j,0.0);
                    }
                 } else if(j == i+1){
                    if(i == doccpi_[0]-2 || i == doccpi_[0]){
                        Rotation->set(0,i,j,-sin(PI));
                    } else{
                        Rotation->set(0,i,j,0.0);
                    }
                } else{
                    Rotation->set(0,i,j,0.0);
                }

            }
        }
        for(int h=1; h < nirrep_; ++h){
            for(int i =0; i < nmopi_[h]; ++i){
                Rotation->set(h,i,i,1.0);
            }
        }
        Rotation->print();

        XCa->gemm(false,false,1.0,Ca_,Rotation,0.0);*/

        FDiaOccA = std::shared_ptr<Vector>(        
        new Vector("Fock matrix occupied diagonal", adoccpi_));
        FDiaVirA = std::shared_ptr<Vector>(        
        new Vector("Fock matrix virtual diagonal", avirtpi_));

//        FockA = Fa();
        Fa_->print();
        Ca_->print();
//        Ca_ = XCa;
//        Ca_->print();
    }
    else{
        throw PSIEXCEPTION("OACCD only implemented for RHF");
    }


}

double Oaccd::compute_energy()
{
    dpdbuf4 T2size; //T2 amplitudes

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

    //Set up DIIS manager, size is set in mp2_energy
    t2DiisManager = new DIISManager(cc_maxdiis, "CCD DIIS T2 amps", DIISManager::LargestError, 
                                    DIISManager::InCore);

    //Start orbital iteration loop here 
    
    int_trans_rhf();
    mp2_energy = mp2_energy_rhf();

    outfile->Printf("\nMP2 energy:  %15.9f \n", mp2_energy);
    outfile->Printf("Total energy:  %15.9f \n", mp2_energy + energy_);
                   
    ccd_energy = ccd_energy_rhf();

    outfile->Printf("\nTotal SCF energy:  %16.10f \n", energy_);
    outfile->Printf("Total MP2 energy:  %16.10f \n", mp2_energy + energy_);
    outfile->Printf("Total CCD energy:  %16.10f \n", ccd_energy + energy_);
                   
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

