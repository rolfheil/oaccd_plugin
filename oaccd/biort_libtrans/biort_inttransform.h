/*
 * @BEGIN LICENSE
 *
 * oaccd by Psi4 Developer, a plugin to:
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

#ifndef BIORTINTTRANS_H
#define BIORTINTTRANS_H

#include "../biort_wfn/biortwfn.h"
#include "psi4/libtrans/integraltransform.h"

using namespace std;

namespace psi{ namespace oaccd {

  /**
   Integral transform to handle biorthogonal transformations.
   Based on IntegralTransform class from Psi4. See this class 
   for additional use instructions
   */

class BiortIntTransform : public IntegralTransform
{
public:
    BiortIntTransform(std::shared_ptr<Wavefunction> wfn,
                          SpaceVec spaces,
                          SharedMatrix lCa,
                          SharedMatrix rCa,
                          SharedMatrix lCb,
                          SharedMatrix rCb,
                          TransformationType transformationType = Restricted,
                          OutputType outputType = DPDOnly,
                          MOOrdering moOrdering = QTOrder,
                          FrozenOrbitals frozenOrbitals = OccAndVir,
                          bool initialize = true);
    
    void initialize();

    void transform_tei(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2,
                          const std::shared_ptr<MOSpace> s3, const std::shared_ptr<MOSpace> s4,
                          HalfTrans ht);

    void transform_tei_first_half(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2);

    void transform_tei_second_half(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2,
                          const std::shared_ptr<MOSpace> s3, const std::shared_ptr<MOSpace> s4);

    void set_orbitals(SharedMatrix C);

    void update_orbitals();   

    virtual ~BiortIntTransform();

protected:

    //Left and right MO coefficient matrices alpha
    std::shared_ptr<Matrix> lCa_;
    std::shared_ptr<Matrix> rCa_;

    //Left and right MO coefficient matrices beta
    std::shared_ptr<Matrix> lCb_;
    std::shared_ptr<Matrix> rCb_;

    // The alpha MO coefficients for all unique spaces needed
    std::map<char, SharedMatrix> laMOCoefficients_;
    // The alpha MO coefficients for all unique spaces needed
    std::map<char, SharedMatrix> raMOCoefficients_;

    // The beta MO coefficients for all unique spaces needed
    std::map<char, SharedMatrix> lbMOCoefficients_;
    // The beta MO coefficients for all unique spaces needed
    std::map<char, SharedMatrix> rbMOCoefficients_;

    void process_eigenvectors();
};

}} // End namespaces

#endif 
