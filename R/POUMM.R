# POUMM.R
# ThreePointUsingSPLITT
# 
# Copyright 2017 Venelin Mitov
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


#' Calculate the POUMM log-likelihood for a given tree, data and model parameters
#' @param x a numerical vector of size N, where N is the number of tips in tree
#' @param tree a phylo object
#' @param x0,alpha,theta,sigma2,sigmae2 parameters of the PMM:
#' \describe{
#' \item{x0}{value at the root of tree excluding whte noise;}
#' \item{alpha}{selection strength parameter of the OU process;}
#' \item{theta}{long-term optimum of the OU process;}
#' \item{sigma2}{unit-time variance increment of the heritable component (OU process);}
#' \item{sigmae2}{variance of the non-heritable component.}
#' }
#' @param ord integer vector: indices of the rows in tree$edge in pruning order. 
#' 
#' Default: \code{reorder(tree, order = "postorder", index.only = TRUE)}.
#' @return the log-likelihood value.
POUMMLogLik <- function(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                        ord = reorder(tree, order = "postorder", index.only = TRUE)) {
  # number of tips in the tree
  N <- length(tree$tip.label)
  # total number of nodes in the tree (tips, internal nodes and root node)
  M <- nrow(tree$edge) + 1L
  # state variables for each node
  a <- b <- c <- rep(0.0, M)
   
  
  for(o in ord) {
    # daughter node
    i <- tree$edge[o, 2]
    # parent node
    j <- tree$edge[o, 1]
    # branch length
    t <- tree$edge.length[o]

    if(i <= N) {
      # initialize a tip node
      a[i] <- -0.5 / sigmae2
      b[i] <- (x[i] - theta) / sigmae2
      c[i] <- -0.5 * ((x[i] - theta) * b[i] + log(2*pi*sigmae2))
    }

    talpha <- t * alpha
    etalpha <- exp(talpha)
    e2talpha <- etalpha * etalpha
    
    if(alpha != 0) {
      u_alpha2t <- alpha / (1 - e2talpha)
    } else {
      u_alpha2t = -0.5 / t
    }
    d <- e2talpha + (a[i] * sigma2) / u_alpha2t
    
    a[j] <- a[j] + a[i] / d
    b[j] <- b[j] + (etalpha * b[i]) / d
    c[j] <- c[j] + -0.5 * log(d) - 0.25 * sigma2 * b[i] * b[i] /
      (u_alpha2t - alpha + a[i] * sigma2) + talpha + c[i]
  }

  # for phylo objects, N+1 denotes the root node
  a[N+1]*(x0 - theta)^2 + b[N+1]*(x0-theta) + c[N+1]
}

#' Calculate the PMM log-likelihood for a given tree, data and model parameters using the RCPP_PMM module
#' @inheritParams POUMMLogLik
#' @param cppObject a previously created object returned by either 
#' \code{\link{New3PointPOUMMCppObject}} or \code{\link{NewAbcPOUMMCppObject}} .
#' @param mode an integer denoting the mode for traversing the tree, i.e. serial vs parallel.
#' 
#' @return the log-likelihood value.
POUMMLogLikCpp <- function(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                         cppObject = New3PointPOUMMCppObject(x, tree),
                         mode = getOption("SPLITT.postorder.mode", 0)) {
  cppObject$TraverseTree(c(x0, alpha, theta, sigma2, sigmae2), mode)
}


#' Create an instance of the RCPP_PMM module for a given tree and trait data
#'
#' @inheritParams POUMMLogLik
#' @return an object to be passed as argument of the \link{POUMMLogLikCpp} function.
#' @seealso \code{\link{PMMLogLikCpp}}
New3PointPOUMMCppObject <- function(x, tree) {
  ThreePointUsingSPLITT__TraversalTaskThreePointPOUMM$new(tree, x[1:length(tree$tip.label)])
}

#' Create an instance of the RCPP_PMM module for a given tree and trait data
#'
#' @inheritParams POUMMLogLik
#' @return an object to be passed as argument of the \link{POUMMLogLikCpp} function.
#' @seealso \code{\link{PMMLogLikCpp}}
NewAbcPOUMMCppObject <- function(x, tree) {
  ThreePointUsingSPLITT__TraversalTaskAbcPOUMM$new(tree, x[1:length(tree$tip.label)])
}
