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
    if(!alreadyPresorted_)
        presort_so_tei();

    // Some
    int nBuckets, thisBucketRows;
    size_t rowsPerBucket, rowsLeft, memFree;

    // Form the Density matrices associated with the C matrices, and allocate the F matrices.
        
    if(Clmat->rowspi() != sopi_)
        throw PSIEXCEPTION("Row dimension of C_l matrix is not equal to SOs per irrep in LibTrans::compute_fock_like_matrices()");
    if(Crmat->rowspi() != sopi_)
        throw PSIEXCEPTION("Row dimension of C_r matrix is not equal to SOs per irrep in LibTrans::compute_fock_like_matrices()");
    SharedMatrix Fmat(new Matrix("F matrix", sopi_, sopi_));
    SharedMatrix Dmat(new Matrix("D matrix", sopi_, sopi_));
    Dmat->gemm(false, true, 1.0, Crmat, Clmat, 0.0);

    outfile->Printf("The density matrix from Biort");
    Dmat->print();

    psio_->open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);

    // Grab control of DPD for now, but store the active number to restore it later
    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    dpdbuf4 J;
    global_dpd_->buf4_init(&J, PSIF_SO_PRESORT, 0, DPD_ID("[n,n]"), DPD_ID("[n,n]"),
                           DPD_ID("[n>=n]+"), DPD_ID("[n>=n]+"), 0, "SO Ints (nn|nn)");

    for(int h=0; h < nirreps_; h++) {
        if(J.params->coltot[h] && J.params->rowtot[h]) {
            memFree = static_cast<size_t>(dpd_memfree() - J.params->coltot[h]);
            rowsPerBucket = memFree / (J.params->coltot[h]);
            if(rowsPerBucket > J.params->rowtot[h]) rowsPerBucket = (size_t) J.params->rowtot[h];
            nBuckets = static_cast<int>(std::ceil(static_cast<double>(J.params->rowtot[h])/
                                             static_cast<double>(rowsPerBucket)));
            rowsLeft = static_cast<size_t>(J.params->rowtot[h] % rowsPerBucket);
        }else{
            nBuckets = 0;
            rowsPerBucket = 0;
            rowsLeft = 0;
        }

        if(print_ > 1) {
            outfile->Printf( "\th = %d; memfree         = %lu\n", h, memFree);
            outfile->Printf( "\th = %d; rows_per_bucket = %lu\n", h, rowsPerBucket);
            outfile->Printf( "\th = %d; rows_left       = %lu\n", h, rowsLeft);
            outfile->Printf( "\th = %d; nbuckets        = %d\n", h, nBuckets);

        }

        // We assume that we're always building totally symmetric Fock matrices here.
        int sym = 0;

        for(int n=0; n < nBuckets; n++){
            if(nBuckets == 1)
                thisBucketRows = rowsPerBucket;
            else
                thisBucketRows = (n < nBuckets-1) ? rowsPerBucket : rowsLeft;
            global_dpd_->buf4_mat_irrep_rd_block(&J, h, n*rowsPerBucket, thisBucketRows);
            for(int pq=0; pq < thisBucketRows; pq++) {
                int pabs = J.params->roworb[h][pq][0];
                int qabs = J.params->roworb[h][pq][1];
                int psym = J.params->psym[pabs];
                int qsym = J.params->qsym[qabs];
                int prel = pabs - J.params->poff[psym];
                int qrel = qabs - J.params->qoff[qsym];
                for(int rs = 0; rs < J.params->coltot[h]; ++rs){
                    int rabs = J.params->colorb[h][rs][0];
                    int sabs = J.params->colorb[h][rs][1];
                    int rsym = J.params->rsym[rabs];
                    int ssym = J.params->ssym[sabs];
                    int rrel = rabs - J.params->roff[rsym];
                    int srel = sabs - J.params->soff[ssym];
                    int pqsym = psym ^ qsym;
                    int rssym = rsym ^ ssym;
                    int qrsym = qsym ^ rsym;
                    int pssym = psym ^ ssym;
                    int prsym = psym ^ rsym;
                    int qssym = qsym ^ ssym;
                    double value = J.matrix[h][pq][rs];

                    if(pqsym == rssym && pqsym == sym){
                        Fmat->add(rsym, rrel, srel, Dmat->get(psym, prel, qrel) * value);
                    }

                    if(qrsym == pssym && qrsym == sym){
                        Fmat->add(qsym, qrel, rrel, -0.5*Dmat->get(psym, prel, srel) * value);
                    }
                } /* rs */
            } /* pq */
        }
        global_dpd_->buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
    } /* h */

    Fmat->add(Hcore);

    global_dpd_->buf4_close(&J);

    psio_->close(PSIF_SO_PRESORT, keepDpdSoInts_);

    // Hand DPD control back to the user
    dpd_set_default(currentActiveDPD);

    return Fmat;
}

/**
 * Presort the two-electron integrals into DPD buffers to prepare them for
 * the transformation.  The frozen core operator is built simultaneously.
 * If this action has already been performed, it will just load the frozen
 * core operator from disk and return.
 */

}}//end namespaces
