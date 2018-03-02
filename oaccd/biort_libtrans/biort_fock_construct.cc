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
SharedMatrix BiortIntTransform::compute_biort_fock_matrix(SharedMatrix Hcore, SharedMatrix Clmat, SharedMatrix Crmat)
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

    // Grab control of DPD for now, but store the active number to restore it later
    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    //Form the density matrices
    std::vector<int> sopiv(nirreps_);
    std::vector<int> occpiv(nirreps_);
    std::vector<unsigned long> zerov(nirreps_);

    for(int h = 0; h < nirreps_; h++){
        sopiv[h] = sopi_[h];
        occpiv[h] = clsdpi_[h];  
        zerov[h]= 0;
    }
    Dmat->gemm('N', 'T', sopiv,sopiv,occpiv,1.0,rCa_, sopiv,lCa_,sopiv,0.0,sopiv);
   
    outfile->Printf("The density matrix swarm");
    Dmat->print();

    //Open the DPD file with integrals
    dpdbuf4 J;
    global_dpd_->buf4_init(&J, PSIF_SO_PRESORT, 0, DPD_ID("[n,n]"), DPD_ID("[n,n]"),
                           DPD_ID("[n>=n]+"), DPD_ID("[n>=n]+"), 0, "SO Ints (nn|nn)");

    for(int h = 0; h < nirreps_; h++){
       outfile->Printf("h: %3i, nirreps: %3i \n",h,nirreps_);
       for(int h_delta = 0; h_delta < nirreps_; h_delta++){
           int h_gamma = h^h_delta;
           if(sopiv[h_gamma] &&  sopiv[h_delta]){
               global_dpd_->buf4_mat_irrep_init_block(&J, h, sopiv[h_gamma]);
               for(int delta = 0; delta < sopiv[h_delta]; delta++){

                   outfile->Printf("h: %3i, h_delta: %3i, delta: %3i \n",h,h_delta,delta);
                   global_dpd_->buf4_mat_irrep_rd_block(&J, h, delta*sopiv[h_gamma],sopiv[h_gamma]);
                 
                   double **pFmat = Fmat->pointer(h);  
                   double **pDmat = Dmat->pointer(h_gamma);  
                   C_DGEMV('N', J.params->rowtot[h], sopiv[h_gamma], 1.0, &J.matrix[h][0][0], J.params->rowtot[h], pDmat[delta], 1, 1.0, pFmat[0],1);

                   pFmat = Fmat->pointer(h_delta);  
                   pDmat = Dmat->pointer(h_gamma);  
                   C_DGEMV('N', sopiv[h_delta], J.params->rowtot[h_gamma],-0.5, &J.matrix[h][0][0],J.params->rowtot[h_gamma],pDmat[0],1, 1.0, pFmat[delta],1);
               } 
           }
       }  
   }


    global_dpd_->buf4_close(&J);

    psio_->close(PSIF_SO_PRESORT, keepDpdSoInts_);

    // Hand DPD control back to the user
    dpd_set_default(currentActiveDPD);

    outfile->Printf("Hcore swarm");
    Hcore->print();

    outfile->Printf("The transformed Fock matrix swarm");
    Fmat->print();

    return Fmat;
}

}}//end namespaces
