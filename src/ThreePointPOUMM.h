/*
 *  ThreePointPOUMM.h
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

#ifndef ThreePointPOUMM_H_
#define ThreePointPOUMM_H_

#include "ThreePointUnivariate.h"
#include "NumericTraitData.h"

using namespace SPLITT;

namespace ThreePointUsingSPLITT {
template<class Tree>
class ThreePointPOUMM: public ThreePointUnivariate<Tree> {

public:
  typedef ThreePointPOUMM<Tree> MyType;
  typedef Tree TreeType;
  typedef PostOrderTraversal<MyType> AlgorithmType;
  typedef ThreePointUnivariate<TreeType> BaseType;
  typedef vec ParameterType;
  typedef NumericTraitData<typename TreeType::NodeType> DataType;
  typedef vec StateType;

  // univariate trait vector
  SPLITT::vec x;
  
  // tree height (maximum root-tip distance)
  double T; 
  // h: height (distance from the root) for each node in the tree
  SPLITT::vec h;
  // u: distance from the far-most tip for each node (, i.e. u[i] = T - h[i])
  SPLITT::vec u;
  
  double x0, alpha, theta, sigma2, sigmae2, e2alphaT, sum_u;
  
  ThreePointPOUMM(
    TreeType const& tree, DataType const& input_data):
    BaseType(tree) {

    if(input_data.x_.size() != this->ref_tree_.num_tips()) {
      throw std::invalid_argument("ERR:01201:SPLITT:ThreePointPOUMM.h:ThreePointPOUMM:: The vector x must be the same length as the number of tips.");
    } else {

      uvec ordNodes = this->ref_tree_.OrderNodes(input_data.names_);
      this->x = At(input_data.x_, ordNodes);
      vec X(this->ref_tree_.num_tips());
      this->set_X_and_Y(X, X);

      // A root-to-node distance vector in the order of pruning processing
      h.resize(this->ref_tree_.num_nodes() - 1);
      for(int i = this->ref_tree_.num_nodes() - 2; i >= 0; i--) {
        h[i] = h[this->ref_tree_.FindIdOfParent(i)] + this->ref_tree_.LengthOfBranch(i);
      }

      this->T = *std::max_element(h.begin(), h.begin()+this->ref_tree_.num_tips());
      
      this->u = SPLITT::vec(this->ref_tree_.num_tips());
      for(int i = 0; i < this->ref_tree_.num_tips(); i++) {
        u[i] = T - h[i];
      }
      sum_u = 0;
      for(auto uu : u) sum_u += uu;
    }
  }

  void SetParameter(ParameterType const& par) {
    if(par.size() != 5) {
      throw std::invalid_argument(
          "ERR:01211:SPLITT:ThreePointPOUMM.h:SetParameter:: The par vector should be of length 5 with \
      elements corresponding to x0, alpha, theta, sigma2 and sigmae2.");
    }
    if(par[1] < 0 || par[3] < 0 || par[4] < 0) {
      throw std::logic_error("ERR:01212:SPLITT:ThreePointPOUMM.h:SetParameter:: The parameters alpha, sigma2 and sigmae2 should be non-negative.");
    }
    this->x0 = par[0];
    this->alpha = par[1];
    this->theta = par[2];
    this->sigma2 = par[3];
    this->sigmae2 = par[4];
    this->e2alphaT = exp(-2*alpha*T);
  }

  inline void InitNode(uint i) {
    ThreePointUnivariate<TreeType>::InitNode(i);
    if(i < this->ref_tree_.num_nodes() - 1) {
      // if an internal node or a tip, transform the branch length leading to this tip
      uint iParent = this->ref_tree_.FindIdOfParent(i);
      
      double ealphahi = exp(alpha*h[i]);

      this->tTransf[i] = sigma2/(2*alpha) *
        (e2alphaT*(ealphahi*ealphahi - exp((2*alpha)*h[iParent])));

      if(i < this->ref_tree_.num_tips()) {
        double mu = theta + (x0 - theta) / ealphahi;
        double ealphaui = exp(-alpha*u[i]);
        this->X[i] = this->Y[i] = (x[i] - mu)*ealphaui;
        this->tTransf[i] += sigmae2 * ealphaui*ealphaui;
      }
    }
  }

  // StateType StateAtRoot() const {
  //   vec res = BaseType::StateAtRoot();
  //   res.push_back(alpha*sum_u);
  //   return res;
  // }
  // 
  inline StateType StateAtRoot() const {
    vec res(1);
    double lnDetVRoot = 2*alpha*sum_u + this->lnDetV[this->ref_tree_.num_nodes() - 1];
    double QRoot = this->Q[this->ref_tree_.num_nodes() - 1];
    res[0] = -0.5*(this->ref_tree_.num_tips() * log(2*G_PI)+ lnDetVRoot + QRoot);
    return res;
  }
};

}
#endif // ThreePointPOUMM_H_
