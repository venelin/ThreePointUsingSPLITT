/*
 *  ThreePointPMM.h
 *  SPLITT
 *
 * Copyright 2017 Venelin Mitov
 *
 * This file is part of SPLITT: a generic C++ library for Serial and Parallel
 * Lineage Traversal of Trees.
 *
 * SPLITT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * SPLITT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SPLITT.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * @author Venelin Mitov
 */

#ifndef ThreePointPMM_H_
#define ThreePointPMM_H_

#include "ThreePointUnivariate.h"
#include "NumericTraitData.h"

using namespace SPLITT;

namespace ThreePointUsingSPLITT {

template<class Tree>
class ThreePointPMM: public ThreePointUnivariate<Tree> {

public:
  typedef ThreePointPMM<Tree> MyType;
  typedef Tree TreeType;
  typedef PostOrderTraversal<MyType> AlgorithmType;
  typedef ThreePointUnivariate<TreeType> BaseType;
  typedef vec ParameterType;
  typedef NumericTraitData<typename TreeType::NodeType> DataType;
  typedef vec StateType;

  // univariate trait vector
  SPLITT::vec x;
  double x0, sigma2, sigmae2;


  ThreePointPMM(
    TreeType const& tree, DataType const& input_data):
    BaseType(tree) {

    if(input_data.x_.size() != this->ref_tree_.num_tips()) {
      throw std::invalid_argument("ERR:01201:SPLITT:ThreePointPMM.h:ThreePointPMM:: The vector x must be the same length as the number of tips.");
    } else {

      uvec ordNodes = this->ref_tree_.OrderNodes(input_data.names_);
      this->x = At(input_data.x_, ordNodes);
      vec X(this->ref_tree_.num_tips());
      this->set_X_and_Y(X, X);
    }
  }

  void SetParameter(ParameterType const& par) {
    if(par.size() != 3) {
      throw std::invalid_argument(
          "ERR:01211:SPLITT:ThreePointPMM.h:SetParameter:: The par vector should be of length 3 with \
      elements corresponding to x0, sigma and sigmae.");
    }
    if(par[1] <= 0 || par[2] <= 0 ) {
      throw std::logic_error("ERR:01212:SPLITT:ThreePointPMM.h:SetParameter:: The parameters sigma and sigmae should be positive.");
    }
    this->x0 = par[0];
    this->sigma2 = par[1];
    this->sigmae2 = par[2];
  }

  inline void InitNode(uint i) {
    ThreePointUnivariate<TreeType>::InitNode(i);
    
    if(i < this->ref_tree_.num_nodes() - 1) {
      // if an internal node or a tip, transform the branch length leading to this tip
      // The call to the parent's class InitNode method has set tTransf[i] to the branch-length.
      this->tTransf[i] = sigma2 * this->tTransf[i];
      if(i < this->ref_tree_.num_tips()) {
        // if an a tip, transform the branch length leading to this tip
        this->X[i] = this->Y[i] = x[i] - x0;
        this->tTransf[i] += sigmae2;
      }  
    }
  }
  
  inline StateType StateAtRoot() const {
    vec res(1);
    double lnDetVRoot = this->lnDetV[this->ref_tree_.num_nodes() - 1];
    double QRoot = this->Q[this->ref_tree_.num_nodes() - 1];
    res[0] = -0.5*(this->ref_tree_.num_tips() * log(2*G_PI)+lnDetVRoot+QRoot);
    return res;
  }
};

}
#endif // ThreePointPMM_H_