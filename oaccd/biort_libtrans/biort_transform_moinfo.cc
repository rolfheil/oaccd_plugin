/*::
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

#include "biort_inttransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/psifiles.h"

#include <math.h>
#include <ctype.h>
#include <stdio.h>

#define EXTERN
#include "psi4/libdpd/dpd.h"

namespace psi{ namespace oaccd {

void BiortIntTransform::process_eigenvectors()
{
    std::vector<std::shared_ptr<MOSpace> >::const_iterator space;

    if(print_ > 4){
        lCa_->print();
        rCa_->print();
        lCb_->print();
        rCb_->print();
    }

    // N.B. The frozen orbitals have been zeroed, if appropriate
    Dimension focc = frzcpi_;
    Dimension aocc = nalphapi_ - frzcpi_;
    Dimension bocc = nbetapi_ - frzcpi_;
    Dimension avir = mopi_ - nalphapi_ - frzvpi_;
    Dimension bvir = mopi_ - nbetapi_ - frzvpi_;
    Dimension aall = mopi_ - frzcpi_ - frzvpi_;
    Dimension fvir = frzvpi_;
    Dimension ball = mopi_ - frzcpi_ - frzvpi_;
    Dimension zero = Dimension(nirreps_);

    for(space = uniqueSpaces_.begin(); space != uniqueSpaces_.end(); ++space){
        std::shared_ptr<MOSpace> moSpace = *space;
        SharedMatrix lCa,rCa,lCb,rCb;
        if(moSpace->label() == MOSPACE_FZC){
            // This is the frozen occupied space
            lCa = lCa_->get_block({zero,sopi_},{zero,focc});
            rCa = rCa_->get_block({zero,sopi_},{zero,focc});
            lCa->set_name("Left alpha frozen occupied orbitals");
            rCa->set_name("Right alpha frozen occupied orbitals");
            if(transformationType_ != Restricted){
                lCb = lCb_->get_block({zero,sopi_},{zero,focc});
                rCb = rCb_->get_block({zero,sopi_},{zero,focc});
                lCb->set_name("Left beta frozen occupied orbitals");
                rCb->set_name("Right beta frozen occupied orbitals");
            }
        }else if(moSpace->label() == MOSPACE_OCC){
            // This is the occupied space
            lCa = lCa_->get_block({zero,sopi_},{focc,focc + aocc});
            rCa = rCa_->get_block({zero,sopi_},{focc,focc + aocc});
            lCa->set_name("Left alpha occupied orbitals");
            rCa->set_name("Right alpha occupied orbitals");
            if(transformationType_ != Restricted){
                lCb = lCb_->get_block({zero,sopi_},{focc,focc + bocc});
                rCb = rCb_->get_block({zero,sopi_},{focc,focc + bocc});
                lCb->set_name("Left beta occupied orbitals");
                rCb->set_name("Right beta occupied orbitals");
            }
        }else if(moSpace->label() == MOSPACE_ALL){
            // This is the full space, sans frozen orbitals
            lCa = lCa_->get_block({zero,sopi_},{focc,focc + aall});
            rCa = rCa_->get_block({zero,sopi_},{focc,focc + aall});
            lCa->set_name("All left alpha orbitals");
            rCa->set_name("All right alpha orbitals");
            if(transformationType_ != Restricted){
                lCb = lCb_->get_block({zero,sopi_},{focc,focc + ball});
                rCb = rCb_->get_block({zero,sopi_},{focc,focc + ball});
                lCb->set_name("All left beta orbitals");
                rCb->set_name("All right beta orbitals");
            }
        }else if(moSpace->label() == MOSPACE_VIR){
            // This is the virtual space
            if(transformationType_ == Restricted){
                // Take the true virtual orbitals, and then append the SOCC orbitals
                std::vector<SharedMatrix> virandsoc;
                virandsoc.push_back(lCa_->get_block({zero,sopi_},{nalphapi_,nalphapi_ + avir}));
                virandsoc.push_back(rCa_->get_block({zero,sopi_},{nalphapi_,nalphapi_ + avir}));
                virandsoc.push_back(lCa_->get_block({zero,sopi_},{clsdpi_,clsdpi_ + openpi_}));
                virandsoc.push_back(rCa_->get_block({zero,sopi_},{clsdpi_,clsdpi_ + openpi_}));
                lCa = Matrix::horzcat(virandsoc);
                rCa = Matrix::horzcat(virandsoc);
                lCa->set_name("Left alpha virtual orbitals");
                rCa->set_name("Right alpha virtual orbitals");
            }else{
                lCa = lCa_->get_block({zero,sopi_},{nalphapi_,nalphapi_ + avir});
                rCa = rCa_->get_block({zero,sopi_},{nalphapi_,nalphapi_ + avir});
                lCa->set_name("Left alpha virtual orbitals");
                rCa->set_name("Right alpha virtual orbitals");
                lCb = lCb_->get_block({zero,sopi_},{nbetapi_,nbetapi_ + bvir});
                rCb = rCb_->get_block({zero,sopi_},{nbetapi_,nbetapi_ + bvir});
                lCb->set_name("Left beta virtual orbitals");
                rCb->set_name("Right beta virtual orbitals");
            }
        }else if(moSpace->label() == MOSPACE_FZV){
            // This is the frozen virtual space
            lCa = lCa_->get_block({zero,sopi_},{mopi_ - frzvpi_,mopi_});
            rCa = rCa_->get_block({zero,sopi_},{mopi_ - frzvpi_,mopi_});
            lCa->set_name("Left alpha frozen virtual orbitals");
            rCa->set_name("Right alpha frozen virtual orbitals");
            if(transformationType_ != Restricted){
                lCb = lCb_->get_block({zero,sopi_},{mopi_ - frzvpi_,mopi_});
                rCb = rCb_->get_block({zero,sopi_},{mopi_ - frzvpi_,mopi_});
                lCb->set_name("Left beta frozen virtual orbitals");
                rCb->set_name("Right beta frozen virtual orbitals");
            }
        }else if(moSpace->label() == MOSPACE_NIL){
            // Do nothing!
        }else{
            char label = moSpace->label();
            if(aOrbsPI_.count(label) == 0){
                std::string err("Space ");
                err += label;
                err += " has not been initialized properly.";
                throw PSIEXCEPTION(err);
            }

            const std::vector<int> & aorbs = moSpace->aOrbs();
            const std::vector<int> & borbs = moSpace->bOrbs();
            int nAOrbs = aorbs.size();
            int nBOrbs = borbs.size();
            std::string name("Alpha orbitals for space ");
            name += label;
            lCa = SharedMatrix(new Matrix(name, nirreps_, sopi_, aOrbsPI_[label]));
            rCa = SharedMatrix(new Matrix(name, nirreps_, sopi_, aOrbsPI_[label]));
            int mo_offsets[8];
            mo_offsets[0] = 0;
            for(int h = 1; h < nirreps_; ++h)
                mo_offsets[h] = mo_offsets[h-1] + mopi_[h-1];

            for(int h = 0; h < nirreps_; ++h){
                int count = 0;
                for(int n = 0; n < nAOrbs; ++n){
                    int orb = aorbs[n];
                    if(mosym_[orb] == h){
                        for(int so = 0; so < sopi_[h]; ++so){
                            lCa->set(h, so, count, lCa_->get(h, so, orb - mo_offsets[h]));
                            rCa->set(h, so, count, rCa_->get(h, so, orb - mo_offsets[h]));
                        }
                        count++;
                    }
                }
            }
            if(transformationType_ != Restricted){
                name = "Beta orbitals for space " + std::string(1, label);
                lCb = SharedMatrix(new Matrix(name, nirreps_, sopi_, bOrbsPI_[label]));
                rCb = SharedMatrix(new Matrix(name, nirreps_, sopi_, bOrbsPI_[label]));
                for(int h = 0; h < nirreps_; ++h){
                    int count = 0;
                    for(int n = 0; n < nBOrbs; ++n){
                        int orb = borbs[n];
                        if(mosym_[orb] == h){
                            for(int so = 0; so < sopi_[h]; ++so){
                                lCb->set(h, so, count, lCb_->get(h, so, orb - mo_offsets[h]));
                                rCb->set(h, so, count, rCb_->get(h, so, orb - mo_offsets[h]));
                            }
                            count++;
                        }
                    }
                }
            }
        }

        if(transformationType_ == Restricted) lCb = lCa;
        if(transformationType_ == Restricted) rCb = rCa;

        laMOCoefficients_[moSpace->label()] = lCa;
        raMOCoefficients_[moSpace->label()] = rCa;
        lbMOCoefficients_[moSpace->label()] = lCb;
        rbMOCoefficients_[moSpace->label()] = rCb;

        if(print_ > 5){
            outfile->Printf( "Orbitals for space %c:-\n",moSpace->label());
            lCa->print();
            rCa->print();
            if (transformationType_ != Restricted)
                lCb->print();
                rCb->print();
            outfile->Printf( "\n\n");
        }
    }// End loop over spaces

}
}}//End namespaces
