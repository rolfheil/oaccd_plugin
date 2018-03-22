/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "biort_inttransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libmints/matrix.h"
#include "psi4/psifiles.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libtrans/mospace.h"

#include <cmath>

#define EXTERN
#include "psi4/libdpd/dpd.h"

using namespace psi;

namespace psi{ namespace oaccd{
/**
 * @brief - Constructs the Fock matrix in SO basis. Uses the density matrix constructed from the two 
 *          coefficient matrices. This makes it possible to construct the Fock matrix in a nonorthogonal 
 *          basis. Assumes that you can hold n^3 integrals in memory
 * @param Hcore - the SO basis core Hamiltonian contribution to the Fock matrix.
 * @param Clmat - Left C coefficient matrix
 * @param Crmat - Right C coefficient matrix
 * @param tFock_energy - returns the new Fock matrix 
 * @return Transformed Fock matrix in SO basis
 */
SharedMatrix BiortIntTransform::compute_biort_fock_matrix(SharedMatrix Hcore, SharedMatrix Clmat, SharedMatrix Crmat, double &tFock_energy)
{
    // This function is supposed to only be called after an initial presort, but we'll check to make sure.
    
    psio_->open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);

    if(!alreadyPresorted_)
        presort_so_tei();

    dpdbuf4 I;
   
    // Form the Density matrices associated with the C matrices, and allocate the F matrices.
    if(Clmat->rowspi() != sopi_)
        throw PSIEXCEPTION("Row dimension of C_l matrix is not equal to SOs per irrep in BiortIntTransform::compute_biort_fock_matrix()");
    if(Crmat->rowspi() != sopi_)
        throw PSIEXCEPTION("Row dimension of C_r matrix is not equal to SOs per irrep in BiortIntTransform::compute_biort_fock_matrix()");
    //SharedMatrix Fmat;
    SharedMatrix Dmat(new Matrix("D matrix", sopi_, sopi_));
    SharedMatrix Fmat(new Matrix("F matrix", sopi_, sopi_));
    Fmat->set_name("Transformed Fock matrix");
    SharedMatrix temp(new Matrix("temp mat", sopi_, sopi_));

    // Grab control of DPD for now, but store the active number to restore it later
    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    //Set up some std::vectors for advanced gemm call
    std::vector<int> sopiv(nirreps_);
    std::vector<int> occpiv(nirreps_);

    for(int h = 0; h < nirreps_; h++){
        sopiv[h] = sopi_[h];
        occpiv[h] = clsdpi_[h];  
    }

    //Construct the density matrix
    Dmat->gemm('N', 'T', sopiv,sopiv,occpiv,1.0,Crmat, sopiv,Clmat,sopiv,0.0,sopiv);
   
    //Open the DPD file with integrals
    dpdbuf4 J;
    global_dpd_->buf4_init(&J, PSIF_SO_PRESORT, 0, DPD_ID("[n,n]"), DPD_ID("[n,n]"),
                           DPD_ID("[n>=n]+"), DPD_ID("[n>=n]+"), 0, "SO Ints (nn|nn)");

    double **pFmat;
    double **pDmat;
    int delta_off;
    int alpha_off;
    int alpha_off2;
    int h_gamma;

    //Calculate F_alpha,beta = sum_gamma,delta D_gamma_delta*(2(alpha,beta|gamma,delta) - (alpha,gamma|beta,delta))
    //Loop over integral symmetry
    for(int h = 0; h < nirreps_; h++){  
        delta_off = 0;

        //Loop over last index symmetry in integrals
        for(int h_delta = 0; h_delta < nirreps_; h_delta++){
            h_gamma = h^h_delta;

            if(!sopi_[h_delta] || !sopi_[h_gamma]) continue;
           
            //Allocate space for integrals and set pointer to F matrix
            pFmat = Fmat->pointer(h_delta);  
            global_dpd_->buf4_mat_irrep_init_block(&J, h, sopiv[h_gamma]);

            //Calculcate offset for exchange
            alpha_off = 0;
            for(int h_beta=0; h_beta < h_gamma; h_beta++){
                for(int h_alpha=0; h_alpha < nirreps_; h_alpha++){
                    if((h_alpha^h_beta) == h){
                        alpha_off = alpha_off + sopi_[h_alpha]*sopi_[h_beta];
                    }
                }
            }
           
            //Loop over final index in integrals
            for(int delta = 0; delta < sopiv[h_delta]; delta++){
 
                //Read in all alpha, beta, and gamma for delta
                global_dpd_->buf4_mat_irrep_rd_block(&J, h, delta_off, sopiv[h_gamma]);
                delta_off = delta_off + sopi_[h_gamma];

                //Calculate Coulomb part
                alpha_off2 = 0;
                if(h == 0){
                    for(int h_alpha = 0; h_alpha < nirreps_; h_alpha++ ){
                    if(!sopi_[h_alpha]) continue;
                        pDmat = Dmat->pointer(h_alpha);  
                        C_DGEMV('N', sopi_[h_delta], sopi_[h_alpha]*sopi_[h_alpha], 2.0, &J.matrix[h][0][alpha_off2], J.params->rowtot[h],
                                pDmat[0], 1, 1.0,&pFmat[0][delta],sopi_[h_delta]);
                        alpha_off2 = alpha_off2 + sopi_[h_alpha]*sopi_[h_alpha];
                    }//end h_alpha for-loop     
                }//end h==0 if statement

                //Calculate exchange part
                pDmat = Dmat->pointer(h_gamma);  
                for(int gamma = 0; gamma < sopi_[h_gamma]; gamma++){
                    for(int alpha = 0; alpha < sopi_[h_delta]; alpha++){
                        pFmat[delta][alpha] -= C_DDOT(sopi_[h_gamma],&J.matrix[h][gamma][alpha+alpha_off],sopi_[h_delta],pDmat[gamma],1);
                    }//end alpha loop
                }//end gamma loop
                         
               
            }// end delta loop
        } //end h_delta loop
    } //end h loop
    global_dpd_->buf4_close(&J);

    psio_->close(PSIF_SO_PRESORT, keepDpdSoInts_);

    // Hand DPD control back to the user
    dpd_set_default(currentActiveDPD);


    //Calculate transformed Fock energy
    temp->gemm(false,false,1.0,Fmat,Dmat,0.0);
    tFock_energy = temp->trace();
    temp->gemm(false,false,1.0,Hcore,Dmat,0.0);
    tFock_energy += 2*temp->trace();

    Fmat->add(Hcore);

    return Fmat;
}

}}//end namespaces
