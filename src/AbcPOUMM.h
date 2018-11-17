//   Copyright 2017 Venelin Mitov
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
//   limitations under the License.
#ifndef ABC_POUMM_H_
#define ABC_POUMM_H_

#include "./SPLITT.h"
#include "./NumericTraitData.h"
#include <iostream>

namespace ThreePointUsingSPLITT {

using namespace SPLITT;

template<class Tree>
class AbcPOUMM: public TraversalSpecification<Tree> {

public:
  typedef AbcPOUMM<Tree> MyType;
  typedef TraversalSpecification<Tree> BaseType;
  typedef Tree TreeType;
  typedef PostOrderTraversal<MyType> AlgorithmType;
  typedef vec ParameterType;
  typedef NumericTraitData<typename TreeType::NodeType> DataType;
  typedef vec StateType;

  double x0, alpha, theta, sigma2, sigmae2;
  vec x;
  vec a, b, c;

  AbcPOUMM(TreeType const& tree, DataType const& input_data):
    BaseType(tree) {

    if(input_data.x_.size() != this->ref_tree_.num_tips() ) {
      std::ostringstream oss;
      oss<<"The vector x must be the same length as the number of tips ("<<
        this->ref_tree_.num_tips()<<"), but were"<<input_data.x_.size()<<".";
      throw std::invalid_argument(oss.str());
    } else {

      uvec ordNodes = this->ref_tree_.OrderNodes(input_data.names_);
      this->x = At(input_data.x_, ordNodes);
      this->a = vec(this->ref_tree_.num_nodes());
      this->b = vec(this->ref_tree_.num_nodes());
      this->c = vec(this->ref_tree_.num_nodes());
    }
  };

  StateType StateAtRoot() const {
    vec res(1);
    double x0_theta = this->x0-this->theta;
    res[0] = a[this->ref_tree_.num_nodes() - 1] * (x0_theta) * (x0_theta) + 
      b[this->ref_tree_.num_nodes() - 1] * x0_theta + 
      c[this->ref_tree_.num_nodes() - 1];
    return res;
  };

  void SetParameter(ParameterType const& par) {
    if(par.size() != 5) {
      throw std::invalid_argument(
      "The par vector should be of length 5 with \
      elements corresponding to x0, alpha, theta, sigma and sigmae.");
    }
    if(par[1] < 0 || par[3] < 0 || par[4] < 0) {
      throw std::logic_error("The parameters alpha, sigma and sigmae should be non-negative.");
    }
    this->x0 = par[0];
    this->alpha = par[1];
    this->theta = par[2];
    this->sigma2 = par[3];
    this->sigmae2 = par[4];
  }

  inline void InitNode(uint i) {
    if(i < this->ref_tree_.num_tips()) {
      double z1 = x[i] - theta;
      a[i] = -0.5 / sigmae2;
      b[i] = z1 / sigmae2;
      c[i] = -0.5 * (M_LN_2PI  + z1 * b[i] + log(sigmae2));
    } else {
      a[i] = b[i] = c[i] = 0;
    }
  }

  inline void VisitNode(uint i) {

    double t = this->ref_tree_.LengthOfBranch(i);
    double talpha = t * alpha;
    double etalpha = exp(talpha);
    double e2talpha = etalpha * etalpha;
    double fe2talpha;
    if(alpha != 0) {
      fe2talpha = alpha / (1 - e2talpha);
    } else {
      fe2talpha = -0.5 / t;
    }
    double gutalphasigma2 = e2talpha + (a[i] * sigma2) / fe2talpha;

    c[i] = -0.5 * log(gutalphasigma2) - 0.25 * sigma2 * b[i] * b[i] /
      (fe2talpha - alpha + a[i] * sigma2) + talpha + c[i];
    b[i] = (etalpha * b[i]) / gutalphasigma2;
    a[i] /= gutalphasigma2;
  }

  inline void PruneNode(uint i, uint i_parent) {
    a[i_parent] += a[i];
    b[i_parent] += b[i];
    c[i_parent] += c[i];
  }

};
}
#endif //ABC_POUMM_H_
