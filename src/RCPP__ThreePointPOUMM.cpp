/**
  *  RCPP__ThreePointPOUMM.cpp
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

#include <Rcpp.h>
#include "./ThreePointPOUMM.h"
    
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]

using namespace SPLITT;
using namespace ThreePointUsingSPLITT;

typedef TraversalTask<
  ThreePointPOUMM<OrderedTree<uint, double>> > TraversalTaskThreePointPOUMM;



TraversalTaskThreePointPOUMM* CreateTraversalTaskThreePointPOUMM(
    Rcpp::List const& tree, vec const& values) {
  
  Rcpp::IntegerMatrix branches = tree["edge"];
  uvec parents(branches.column(0).begin(), branches.column(0).end());
  uvec daughters(branches.column(1).begin(), branches.column(1).end());
  vec t = Rcpp::as<vec>(tree["edge.length"]);
  uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  uvec tip_names = Seq(uint(1), num_tips);
  
  typename TraversalTaskThreePointPOUMM::DataType data(tip_names, values);
  
  return new TraversalTaskThreePointPOUMM(parents, daughters, t, data);
}


// This will enable returning a copy of the `TraversalAlgorithm`-object stored in
// a `TraversalTaskThreePointPOUMM` object to a R. This will be used in the MiniBenchmark
// R-function to check things like the OpenMP version used during compilation and
// the number of OpenMP threads at runtime. 
RCPP_EXPOSED_CLASS_NODECL(TraversalTaskThreePointPOUMM::AlgorithmType)
  
RCPP_MODULE(ThreePointUsingSPLITT__TraversalTaskThreePointPOUMM) {
  
  // Expose the properties VersionOPENMP and NumOmpThreads from the base 
  // TraversalAlgorithm class
  Rcpp::class_<TraversalTaskThreePointPOUMM::AlgorithmType::ParentType> (
      "ThreePointUsingSPLITT__ThreePointPOUMM__TraversalAlgorithm"
    )
  .property( "VersionOPENMP",
             &TraversalTaskThreePointPOUMM::AlgorithmType::ParentType::VersionOPENMP )
  .property( "NumOmpThreads",
             &TraversalTaskThreePointPOUMM::AlgorithmType::ParentType::NumOmpThreads )
  ;

  // Expose the TraversalTaskThreePointPOUMM::AlgorithmType specifying that it derives 
  // from the base TraversalAlgorithm class
  Rcpp::class_<TraversalTaskThreePointPOUMM::AlgorithmType> (
      "ThreePointUsingSPLITT__ThreePointPOUMM__AlgorithmType"
    )
  .derives<TraversalTaskThreePointPOUMM::AlgorithmType::ParentType>(
      "ThreePointUsingSPLITT__ThreePointPOUMM__TraversalAlgorithm"
    )
  ;
  
  // Finally, expose the TraversalTaskThreePointPOUMM class - this is the main class in 
  // the module, which will be instantiated from R using the factory function
  // we've just written.
  Rcpp::class_<TraversalTaskThreePointPOUMM>( "ThreePointUsingSPLITT__TraversalTaskThreePointPOUMM" )
  // The <argument-type-list> MUST MATCH the arguments of the factory function 
  // defined above.
  .factory<Rcpp::List const&, vec const&>( &CreateTraversalTaskThreePointPOUMM )
  // Expose the method that we will use to execute the TraversalTask
  .method( "TraverseTree", &TraversalTaskThreePointPOUMM::TraverseTree )
  // Expose the algorithm property
  .property( "algorithm", &TraversalTaskThreePointPOUMM::algorithm )
  ;
}

