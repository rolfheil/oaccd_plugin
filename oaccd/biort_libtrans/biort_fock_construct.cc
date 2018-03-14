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
 * @brief Computes Fock matrices, frozen core operators and other Fock-like quantities.  This shouldn't
 *        be needed because those quantities are computed during the SO integral presort.  However, in
 *        Brueckner CC and other orbital optimized methods, the Fock matrices need an update so this
 *        function may be called to extract them from the DPD buffers, allowing the IWL buffers to be
 *        nuked after the first integral transformation..
 * @param Hcore - the SO basis core Hamiltonian contribution to the Fock matrix.
 * @param Cmats - a list of the subset of C matrices for each Fock-like quantity, e.g. the nso*nfzc
 *                subset for the frozen core operator.
 * @return a list of the Fock matrices associated with each incoming C matrix.
 */
SharedMatrix BiortIntTransform::compute_biort_fock_matrix(SharedMatrix Hcore, SharedMatrix Clmat, SharedMatrix Crmat, double &tFock_energy)
{
    // This function is supposed to only be called after an initial presort, but we'll check to make sure.
    
    psio_->open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);

    if(!alreadyPresorted_)
        presort_so_tei();

    // Some
    int nBuckets, thisBucketRows;
    size_t rowsPerBucket, rowsLeft, memFree;
    dpdbuf4 I;
   
    // tFa_ = std::shared_ptr<Matrix>(new Matrix(Fa_));:
    // tFb_ = std::shared_ptr<Matrix>(new Matrix(Fb_));

    // Form the Density matrices associated with the C matrices, and allocate the F matrices.
        
    if(Clmat->rowspi() != sopi_)
        throw PSIEXCEPTION("Row dimension of C_l matrix is not equal to SOs per irrep in LibTrans::compute_fock_like_matrices()");
    if(Crmat->rowspi() != sopi_)
        throw PSIEXCEPTION("Row dimension of C_r matrix is not equal to SOs per irrep in LibTrans::compute_fock_like_matrices()");
    //SharedMatrix Fmat;
    // Fmat = std::shared_ptr<Matrix>(new Matrix(Hcore));
    //  Fmat->set_name("Transformed Fock matrix");
    SharedMatrix Dmat(new Matrix("D matrix", sopi_, sopi_));
    SharedMatrix Fmat(new Matrix("F matrix", sopi_, sopi_));
    SharedMatrix temp(new Matrix("temp mat", sopi_, sopi_));

    // Grab control of DPD for now, but store the active number to restore it later
    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    //Form the density matrices
    std::vector<int> sopiv(nirreps_);
    std::vector<int> occpiv(nirreps_);

    for(int h = 0; h < nirreps_; h++){
        sopiv[h] = sopi_[h];
        occpiv[h] = clsdpi_[h];  
    }
    Crmat->print();
    Clmat->print();
    Dmat->print();
    for(int h = 0; h < nirreps_; h++){
        outfile->Printf("h: %3i,  sopiv: %3i \n",h,sopiv[h]);
        outfile->Printf("h: %3i, occpiv: %3i \n",h,occpiv[h]);
    }
    Dmat->gemm('N', 'T', sopiv,sopiv,occpiv,1.0,Crmat, sopiv,Clmat,sopiv,0.0,sopiv);
   
    outfile->Printf("The density matrix swarm");
    Dmat->print();
/*
    double** pCrmat;
    double** pClmat; 
    double** pDmat;
    
    for(int h = 0; h < nirreps_; h++){
        if(sopi_[h]!= 0 && clsdpi_[h] != 0){
            pCrmat = Crmat->pointer(h);  
            pClmat = Clmat->pointer(h); 
            pDmat = Dmat->pointer(h); 
            C_DGEMM('N', 'T', sopi_[h],sopi_[h],clsdpi_[h],1.0,pCrmat[0], sopi_[h],pClmat[0],sopi_[h],0.0, pDmat[0], sopi_[h]);
        }
    }*/
    outfile->Printf("The density matrix swarms again");
    Dmat->print();

    //Open the DPD file with integrals
    dpdbuf4 J;
    global_dpd_->buf4_init(&J, PSIF_SO_PRESORT, 0, DPD_ID("[n,n]"), DPD_ID("[n,n]"),
                           DPD_ID("[n>=n]+"), DPD_ID("[n>=n]+"), 0, "SO Ints (nn|nn)");

    int h = 0;
    double **pFmat;
    double **pDmat;
    int delta_off;
    int alpha_off;
    int h_gamma;
    outfile->Printf("h: %3i, nirreps: %3i \n",h,nirreps_);

    for(int h = 0; h < nirreps_; h++){
        for(int h_delta = 0; h_delta < nirreps_; h_delta++){
            h_gamma = h^h_delta;

            if(!sopi_[h_delta] || !sopi_[h] || !sopi_[h_gamma]) continue;
           
            delta_off = 0;
            pFmat = Fmat->pointer(h_delta);  
            global_dpd_->buf4_mat_irrep_init_block(&J, h, sopiv[h_gamma]);
           
           for(int delta = 0; delta < sopiv[h_delta]; delta++){
 
               outfile->Printf("h: %3i, h_delta: %3i, delta: %3i, sopiv[h_delta]: %3i, delta_off: %3i \n",h,h_delta,delta, sopiv[h_delta], delta_off);
               global_dpd_->buf4_mat_irrep_rd_block(&J, h, delta_off, sopiv[h_gamma]);
               delta_off = delta_off + sopi_[h_gamma];
               outfile->Printf("\n");
/*             for(int n=0; n<sopi_[h_delta]; n++){
                    for(int m=0; m<J.params->rowtot[h]; m++){
                            outfile->Printf("%10.6f", J.matrix[h][n][m]);
                    }
                    outfile->Printf("\n");
               }
               outfile->Printf("\n");
               outfile->Printf("\n");*/
               outfile->Printf("Before swarm\n");
               for(int n=0; n<sopi_[h_delta]; n++){
                    for(int m=0; m<sopi_[h_delta]; m++){
                            outfile->Printf("%10.6f", pFmat[n][m]);
                    }
                    outfile->Printf("\n");
               }
               outfile->Printf("\n");
               outfile->Printf("\n");

               alpha_off = 0;
               if(h == 0){
                    for(int h_alpha = 0; h_alpha < nirreps_; h_alpha++ ){
                        if(!sopi_[h_alpha]) continue;
                         pDmat = Dmat->pointer(h_alpha);  
                         C_DGEMV('N', sopi_[h_delta], sopi_[h_alpha]*sopi_[h_alpha], 2.0, &J.matrix[h][0][alpha_off], J.params->rowtot[h],
                                  pDmat[0], 1, 1.0,&pFmat[0][delta],sopi_[h_delta]);
                         alpha_off = alpha_off + sopi_[h_alpha]*sopi_[h_alpha];
                    }//end alpha for-loop 
               }//end if statement 
                         
               pDmat = Dmat->pointer(h_gamma);  
               C_DGEMV('T', sopi_[h_gamma]*sopi_[h_gamma], sopi_[h_delta],-1.0, &J.matrix[h][0][0], sopi_[h_delta],pDmat[0], 1, 1.0, pFmat[delta],1);
               
               outfile->Printf("\n");
               outfile->Printf("After swarm \n");
               for(int n=0; n<sopi_[h_delta]; n++){
                    for(int m=0; m<sopi_[h_delta]; m++){
                            outfile->Printf("%10.6f", pFmat[n][m]);
                    }
                    outfile->Printf("\n");
               }
               outfile->Printf("\n");

           }// end delta loop
       } //end h_delta loop
    } //end h loop
    global_dpd_->buf4_close(&J);

    psio_->close(PSIF_SO_PRESORT, keepDpdSoInts_);

    // Hand DPD control back to the user
    dpd_set_default(currentActiveDPD);

    outfile->Printf("Hcore swarm");
    Hcore->print();

    temp->transform(Clmat,Fmat,Crmat);
    double exchange = 0.0;
    for(int h = 0; h < nirreps_; h++){
        for(int i = 0; i < occpiv[h]; i++){
            outfile->Printf("h: %3i, nirreps: %3i \n",h,nirreps_);
            exchange -= temp->get(h,i,i); 
        }
    }
    
    outfile->Printf("\n the Trace is: %f \n",exchange);

    outfile->Printf("The transformed Fock matrix swarm");
    Fmat->print();

    Fmat->add(Hcore);
    outfile->Printf("The transformed Fock matrix swarm");
    Fmat->print();

    temp->transform(Clmat,Fmat,Crmat);
    temp->print();
    for(int h = 0; h < nirreps_; h++){
        for(int i = 0; i < occpiv[h]; i++){
            outfile->Printf("h: %3i, nirreps: %3i \n",h,nirreps_);
            exchange += 2*temp->get(h,i,i); 
        }
    }
    
    tFock_energy = exchange;

    return Fmat;
}

}}//end namespaces
