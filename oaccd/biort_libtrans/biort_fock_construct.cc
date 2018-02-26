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
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    if(!alreadyPresorted_)
        presort_so_tei();

    // Some
    int nBuckets, thisBucketRows;
    size_t rowsPerBucket, rowsLeft, memFree;
    dpdbuf4 I;
    dpdfile2 D, F;
   
    // tFa_ = std::shared_ptr<Matrix>(new Matrix(Fa_));:
    // tFb_ = std::shared_ptr<Matrix>(new Matrix(Fb_));

    // Form the Density matrices associated with the C matrices, and allocate the F matrices.
        
    if(Clmat->rowspi() != sopi_)
        throw PSIEXCEPTION("Row dimension of C_l matrix is not equal to SOs per irrep in LibTrans::compute_fock_like_matrices()");
    if(Crmat->rowspi() != sopi_)
        throw PSIEXCEPTION("Row dimension of C_r matrix is not equal to SOs per irrep in LibTrans::compute_fock_like_matrices()");
    SharedMatrix Fmat(new Matrix("F matrix", sopi_, sopi_));
    SharedMatrix Dmat(new Matrix("D matrix", sopi_, sopi_));

    // Grab control of DPD for now, but store the active number to restore it later
    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    dpdbuf4 J;
    global_dpd_->buf4_init(&J, PSIF_SO_PRESORT, 0, DPD_ID("[n,n]"), DPD_ID("[n,n]"),
                           DPD_ID("[n>=n]+"), DPD_ID("[n>=n]+"), 0, "SO Ints (nn|nn)");


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

    //Make the density matrices DPD
       
    global_dpd_->file2_init(&D, PSIF_LIBTRANS_DPD, 0, DPD_ID('n'),DPD_ID('n'), "D Transformed");
    global_dpd_->file2_mat_init(&D);

    for(int h = 0; h < nirreps_; ++h){
        for(int alpha = 0; alpha < sopi_[h]; alpha++){
            for(int beta = 0; beta < sopi_[h]; beta++){
                D.matrix[h][alpha][beta] = Dmat->get(h,alpha,beta);
            }
        }
    }

    global_dpd_->file2_mat_wrt(&D);
    global_dpd_->file2_close(&D);


    Fmat->add(Hcore);

    global_dpd_->buf4_close(&J);

    psio_->close(PSIF_SO_PRESORT, keepDpdSoInts_);
    psio_->close(PSIF_LIBTRANS_DPD, keepDpdSoInts_);

    // Hand DPD control back to the user
    dpd_set_default(currentActiveDPD);

    return Fmat;
}

}}//end namespaces
