# PMM.R
# ThreePointUsingSPLITT
# 
# Copyright 2018 Venelin Mitov
# 
# This file is part of SPLITT: a generic C++ library for Serial and Parallel
# Lineage Traversal of Trees.
# 
# ThreePointUsingSPLITT is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
# 
# ThreePointUsingSPLITT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with SPLITT.  If not, see
# <http://www.gnu.org/licenses/>.
# 
# @author Venelin Mitov


#' Calculate the PMM log-likelihood for a given tree, data and model parameters using the Rcpp module
#' @param x a numerical vector of size N, where N is the number of tips in tree
#' @param tree a phylo object
#' @param x0,sigma2,sigmae2 parameters of the PMM:
#' \describe{
#' \item{x0}{value at the root of tree excluding white noise;}
#' \item{sigma2}{unit-time variance increment of the heritable component;}
#' \item{sigmae2}{variance of the non-heritable component.}
#' }
#' @param cppObject a previously created object returned by \code{\link{NewPMMCppObject}}
#' @param mode an integer denoting the mode for traversing the tree, i.e. serial vs parallel.
#' 
#' @return the log-likelihood value.
PMMLogLikCpp <- function(x, tree, x0, sigma2, sigmae2, 
                         cppObject = NewPMMCppObject(x, tree),
                         mode = getOption("SPLITT.postorder.mode", 0)) {
  cppObject$TraverseTree(c(x0, sigma2, sigmae2), mode)
}


#' Create an instance of the Rcpp module for a given tree and trait data
#'
#' @inheritParams POUMMLogLik
#' 
#' @return an object to be passed as argument of the \link{POUMMLogLikCpp} function.
#' @seealso \link{POUMMLogLikCpp}
NewPMMCppObject <- function(x, tree) {
  ThreePointUsingSPLITT__TraversalTaskThreePointPMM$new(tree, x[1:length(tree$tip.label)])
}
