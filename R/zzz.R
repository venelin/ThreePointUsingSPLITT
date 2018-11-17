# zzz.R
# PMMUsingSPLITT
# 
# Copyright 2018 Venelin Mitov
# 
# This file is part of SPLITT: a generic C++ library for Serial and Parallel
# Lineage Traversal of Trees.
# 
# SPLITT is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
# 
# SPLITT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with SPLITT.  If not, see
# <http://www.gnu.org/licenses/>.
# 
# @author Venelin Mitov

## Up until R 2.15.0, the require("methods") is needed but (now)
## triggers an warning from R CMD check
#.onLoad <- function(libname, pkgname){
#    #require("methods")  ## needed with R <= 2.15.0
#    loadRcppModules()
#}


#' Rcpp module for the \code{TraversalTaskAbcPOUMM}-class
#' @name ThreePointUsingSPLITT__TraversalTaskAbcPOUMM
#' @aliases Rcpp_ThreePointUsingSPLITT__TraversalTaskAbcPOUMM-class
NULL

#' \code{TraversalAlgorithm}-type used in \code{AbcPOUMM}
#' @name ThreePointUsingSPLITT__AbcPOUMM__AlgorithmType
#' @aliases Rcpp_ThreePointUsingSPLITT__AbcPOUMM__AlgorithmType-class
NULL

#' Base class for \code{ThreePointUsingSPLITT::AbcPOUMM::AlgorithmType}
#' @name ThreePointUsingSPLITT__AbcPOUMM__TraversalAlgorithm
#' @aliases Rcpp_ThreePointUsingSPLITT__AbcPOUMM__TraversalAlgorithm-class


#' Rcpp module for the \code{TraversalTaskThreePointPOUMM}-class
#' @name ThreePointUsingSPLITT__TraversalTaskThreePointPOUMM
#' @aliases Rcpp_ThreePointUsingSPLITT__TraversalTaskThreePointPOUMM-class
NULL

#' \code{TraversalAlgorithm}-type used in \code{ThreePointPOUMM}
#' @name ThreePointUsingSPLITT__ThreePointPOUMM__AlgorithmType
#' @aliases Rcpp_ThreePointUsingSPLITT__ThreePointPOUMM__AlgorithmType-class
NULL

#' Base class for \code{ThreePointUsingSPLITT::ThreePointPOUMM::AlgorithmType}
#' @name ThreePointUsingSPLITT__ThreePointPOUMM__TraversalAlgorithm
#' @aliases Rcpp_ThreePointUsingSPLITT__ThreePointPOUMM__TraversalAlgorithm-class
NULL


#' Rcpp module for the \code{TraversalTaskThreePointPMM}-class
#' @name ThreePointUsingSPLITT__TraversalTaskThreePointPMM
#' @aliases Rcpp_ThreePointUsingSPLITT__TraversalTaskThreePointPMM-class
NULL

#' \code{TraversalAlgorithm}-type used in \code{ThreePointPMM}
#' @name ThreePointUsingSPLITT__ThreePointPMM__AlgorithmType
#' @aliases Rcpp_ThreePointUsingSPLITT__ThreePointPMM__AlgorithmType-class
NULL

#' Base class for \code{ThreePointUsingSPLITT::ThreePointPMM::AlgorithmType}
#' @name ThreePointUsingSPLITT__ThreePointPMM__TraversalAlgorithm
#' @aliases Rcpp_ThreePointUsingSPLITT__ThreePointPMM__TraversalAlgorithm-class
NULL
# loading the RCPP C++ modules

loadModule( "ThreePointUsingSPLITT__TraversalTaskThreePointPMM", TRUE )
loadModule( "ThreePointUsingSPLITT__TraversalTaskThreePointPOUMM", TRUE )
loadModule( "ThreePointUsingSPLITT__TraversalTaskAbcPOUMM", TRUE )
