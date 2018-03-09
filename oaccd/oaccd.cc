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
#include "math.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libmints/molecule.h"

using namespace std;

namespace psi{ namespace oaccd {

//Will eventually become NOCCD and OACCD routine
//Currently working on the CCD part
//This is similar to the CCSD and OCC stuff in Psi4

Oaccd::Oaccd(SharedWavefunction ref_wfn, Options& options)
    : Biortwfn(options)
{
    // Shallow copy ref_wfn data into this wavefunction
    shallow_copy(ref_wfn);
    reference_wavefunction_ = ref_wfn;
    common_init();
}

Oaccd::~Oaccd()
{
}

void Oaccd::common_init()
{
    std::shared_ptr<Matrix> U_p;
    std::shared_ptr<Matrix> U_m;
    std::shared_ptr<Matrix> kappa;
    double theta = 1.0;
    bool biort= false;
    bool nonsym = false;
    bool t1trans = true;
    bool orthogonal = false;

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

    //Set up some reference energy stuff that we'll need when changing basis
    tF_energy = energy_;
    nuc_energy = reference_wavefunction_->molecule()->nuclear_repulsion_energy(reference_wavefunction_->get_dipole_field_strength());
        

    Fa_->print();
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


        FDiaOccA = std::shared_ptr<Vector>(        
        new Vector("Fock matrix occupied diagonal", adoccpi_));
        FDiaVirA = std::shared_ptr<Vector>(        
        new Vector("Fock matrix virtual diagonal", avirtpi_));

        outfile->Printf("Ca oaccd");
        Ca_->print();
        
        lCa_ = std::shared_ptr<Matrix>(
               new Matrix(Ca_));
        lCa_->set_name("left C matrix");

        rCa_ = std::shared_ptr<Matrix>(
               new Matrix(Ca_));
        rCa_->set_name("right C matrix");

        U_p = std::shared_ptr<Matrix>(
               new Matrix(Ca_));

        kappa = std::shared_ptr<Matrix>(
               new Matrix(Ca_));

        U_p->set(0.0);
        kappa->set(0.0);
        
        for(int h = 0; h < nirrep_; h++){
                for(int i = 0; i < nmopi_[h]; i++){
                        U_p->set(h,i,i,1.0);
                }
        }

        if(biort){
            U_p->set(0,1,1,cosh(theta));
            U_p->set(0,2,2,cosh(theta));
      
            U_p->set(0,5,5,cosh(theta));
            U_p->set(0,6,6,cosh(theta));
      
            U_m = std::shared_ptr<Matrix>(
                   new Matrix(U_p));
 
            U_p->set(0,1,2,-sinh(theta));
            U_p->set(0,2,1,-sinh(theta));
            U_m->set(0,1,2,sinh(theta));
            U_m->set(0,2,1,sinh(theta));
 
            U_p->set(0,5,6,-sinh(theta));
            U_p->set(0,6,5,-sinh(theta));
            U_m->set(0,5,6,sinh(theta));
            U_m->set(0,6,5,sinh(theta));
        }
        else if(nonsym){
            double x = 0.1;

            U_m = std::shared_ptr<Matrix>(
                   new Matrix(U_p));
 
            U_p->set(0,1,2,x);
            U_m->set(0,1,2,-x);
        }
        else if(t1trans){

            U_m = std::shared_ptr<Matrix>(
                   new Matrix(U_p));

/*            kappa->set(0, 3,0, 0.000023187740);
            kappa->set(0, 3,1,-0.000988382087);
            kappa->set(0, 3,2, 0.010371795264);

            kappa->set(2, 1,0,-0.003295502210);*/

            /*
            kappa->set(0, 5,0,-0.000108178);
            kappa->set(0, 5,1,-0.003184533);
            kappa->set(0, 5,2, 0.000000000);
            kappa->set(0, 5,3, 0.007586514);
            kappa->set(0, 5,4, 0.000000000);
            kappa->set(0, 6,0, 0.000000000);
            kappa->set(0, 6,1, 0.000000000);
            kappa->set(0, 6,2,-0.005945744);
            kappa->set(0, 6,3, 0.000000000);
            kappa->set(0, 6,4, 0.000000000);
            kappa->set(0, 7,0, 0.000000000);
            kappa->set(0, 7,1,-0.000000000);
            kappa->set(0, 7,2,-0.001350703);
            kappa->set(0, 7,3, 0.000000000);
            kappa->set(0, 7,4, 0.000000000);
            kappa->set(0, 8,0,-0.000075723);
            kappa->set(0, 8,1, 0.000086889);
            kappa->set(0, 8,2, 0.000000000);
            kappa->set(0, 8,3, 0.001206825);
            kappa->set(0, 8,4, 0.000000000);
            kappa->set(0, 9,0, 0.000000000);
            kappa->set(0, 9,1, 0.000000000);
            kappa->set(0, 9,2, 0.000000000);
            kappa->set(0, 9,3,-0.000000000);
            kappa->set(0, 9,4, 0.000942717);
            kappa->set(0,10,0,-0.000079432);
            kappa->set(0,10,1,-0.003453193);
            kappa->set(0,10,2,-0.000000000);
            kappa->set(0,10,3, 0.003944658);
            kappa->set(0,10,4,-0.000000000);
            kappa->set(0,11,0, 0.000000000);
            kappa->set(0,11,1, 0.000000000);
            kappa->set(0,11,2, 0.008504124);
            kappa->set(0,11,3,-0.000000000);
            kappa->set(0,11,4, 0.000000000);
            kappa->set(0,12,0,-0.000241174);
            kappa->set(0,12,1, 0.002523827);
            kappa->set(0,12,2, 0.000000000);
            kappa->set(0,12,3,-0.002107662);
            kappa->set(0,12,4, 0.000000000);
            */
 
            U_p->axpy(1.0,kappa);
            U_m->axpy(-1.0,kappa);
 
        }
        else if(orthogonal){
            U_p->set(0,1,1,cos(theta));
            U_p->set(0,2,2,cos(theta));
      
            U_m = std::shared_ptr<Matrix>(
                   new Matrix(U_p));
 
            U_p->set(0,1,2,-sin(theta));
            U_p->set(0,2,1,sin(theta));
            U_m->set(0,1,2,sin(theta));
            U_m->set(0,2,1,-sin(theta));
        }
        else{
            U_m = std::shared_ptr<Matrix>(
                   new Matrix(U_p));
        }

        rCb_ = std::shared_ptr<Matrix>(
               new Matrix(Ca_));
        
        rCb_->gemm(false, false, 1.0, rCa_, U_p, 0.0); 
        rCa_->copy(rCb_); 

        rCb_->gemm(false, true, 1.0, lCa_, U_m, 0.0); 
        lCa_->copy(rCb_); 
        
        outfile->Printf("Original Ca_");
        Ca_->print();
        outfile->Printf("lCa_");
        lCa_->print();
        outfile->Printf("rCa_");
        rCa_->print();

        outfile->Printf("swarms \n");

        lCb_ = std::shared_ptr<Matrix>(
               new Matrix(S_));
        outfile->Printf("S overlap matrix");
        lCb_->print();

        lCb_->gemm(false, false, 1.0, S_, lCa_, 0.0);

        rCb_->gemm(true, false, 1.0, lCb_, rCa_, 0.0);
        outfile->Printf("\n Product of S x lCa x rCa \n");
        rCb_->print();

        rCb_->gemm(true, false, 1.0, lCa_, rCa_, 0.0);
        outfile->Printf("\n Product of lCa x rCa \n");
        rCb_->print();

/*        rCb_->gemm('N', 'T', nsopi_[0],nsopi_[0],doccpi_[0],1.0,Ca_,nsopi_[0],Ca_,nsopi_[0],0.0,
                    nsopi_[0],0,0,0);
        outfile->Printf("\n Product of Ca_o x Ca_o \n");
        rCb_->print();
        rCb_->gemm('N', 'T', nsopi_[0],nsopi_[0],doccpi_[0],1.0,rCa_,nsopi_[0],lCa_,nsopi_[0],0.0,
                    nsopi_[0],0,0,0);
        outfile->Printf("\n Product of rCa_o x lCa_o \n");
        rCb_->print();*/

        outfile->Printf("\n Alpha density Da_ \n");
        Da_->print();

        outfile->Printf("\n Core Hamiltonian H_ \n");
        H_->print();

        outfile->Printf("\n MO Fock matrix \n");
        lCb_->transform(Ca_, Fa_, Ca_);
        lCb_->print();

        outfile->Printf("\n U_p matrix \n");
        U_p->print();
        outfile->Printf("\n U_m matrix \n");
        U_m->print();
        
        lCb_ = lCa_;
        rCb_ = rCa_;
        Ca_ ->set(0.0);
    }
    else{
        throw PSIEXCEPTION("OACCD only implemented for RHF");
    }


}

double Oaccd::compute_energy()
{
    //Start by getting the required integrals
    
    //Allocate integrals,must be done after constructor
    const bool initialize=false;

    std::vector<std::shared_ptr<MOSpace>> spaces = {MOSpace::occ, MOSpace::vir};
    ints = new BiortIntTransform(shared_from_this(), spaces,
               lCa(), rCa(), lCb(), rCb(),
               IntegralTransform::TransformationType::Restricted,
               IntegralTransform::OutputType::DPDOnly,
               IntegralTransform::MOOrdering::QTOrder,
               IntegralTransform::FrozenOrbitals::None,
               initialize);

    ints->set_dpd_id(0);    
    ints->set_keep_dpd_so_ints(true);

    outfile->Printf("\n Swarm here and now \n\n");

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

    outfile->Printf("\nTotal SCF energy:  %16.14f \n", energy_);
    outfile->Printf("Trans ref energy:  %16.14f \n", tF_energy);
    outfile->Printf("Total MP2 energy:  %16.14f \n", mp2_energy + tF_energy);
    outfile->Printf("Total CCD energy:  %16.14f \n", ccd_energy + tF_energy);
                   
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

