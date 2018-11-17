
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis build status](https://travis-ci.org/venelin/ThreePointUsingSPLITT.svg?branch=master)](https://travis-ci.org/venelin/ThreePointUsingSPLITT) [![Coverage status](https://codecov.io/gh/venelin/ThreePointUsingSPLITT/branch/master/graph/badge.svg)](https://codecov.io/github/venelin/ThreePointUsingSPLITT?branch=master) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/ThreePointUsingSPLITT?color=blue)](http://cran.r-project.org/web/packages/ThreePointUsingSPLITT)

ThreePointUsingSPLITT
=====================

The 3-point structure algorithm is a linear time algorithm for calculating the determinant, |*V*|, and any quadratic product of the form *a*<sup>*T*</sup>*V*<sup>−1</sup>*b*, where *V* is a *N* × *N* covariance matrix satisfying a 3-point structure, and a is a column vector of size N. This algorithm has been introduced in the article "A linear-time algorithm for Gaussian and non-Gaussian trait evolution models" by Lam si Tung Ho and Cecile Ane published in Systematic Biology in 2014. A matrix V has a 3-point structure if it is symmetric, with nonnegative entries (*V*<sub>*i**j*</sub> &gt; =0) and satisfies the following 3-point condition: for any *i*, *j*, *k* (not necessarily distinct), the two smallest of V\_{ij} , V\_{ik} , and V\_{jk} are equal (Ho and Ane, 2014). This 3-point structure property implies that the matrix *V* can be encoded in a tree of *N* tips, where any element V\_{ij} equals the distance from the root to the most recent common ancestor of the tips *i* and *j* with "most recent common ancestor" meaning the node farthest away from the root of the tree that is an ancestor of both *i* and *j*. This is equaivalent to the covariance structure arising from a univariate branching Wiener stochastic process. The algorithm is generalized to multivaraite such processes (Ho and Ane, 2014). The algorithm does not need to construct the matrix *V* in the memory, hence it is *O*(*N*) complex, both in terms of time and memory. This R-package implements parallel version of this algorithm, based on the 'SPLITT' library for parallel lineage traversal of trees <https://venelin.github.io/SPLITT/index.html>. The header file 'ThreePointUnivariate.h' provides the univariate implementation; the header 'ThreePointMultivariate.h' provides the multivariate implementation. These two headers implement traversal specifications as base C++ classes that can be inherited by specific model classes. Two examples of such model classes are the PMM model implemented in the header 'ThreePointPMM.h' and the POUMM model implemented in the header file 'ThreePointUnivariate.h'. The purpose of these two classes is to provide fast parallel likelihood evaluation for the two models given a phylogenetic tree and univariate trait data measured at its tips. Both of the classes can be invoked from R via Rcpp modules based on the package Rcpp <https://CRAN.R-project.org/package=Rcpp>.

The goal of ThreePointUsingSPLITT is to provide an example of how to use the SPLITT C++ library in an R-package. The package implements parallelized log-likelihood calculation of the univariate phylogenetic mixed model (PMM) and the univariate phylogenetic Ornstein-Uhlenbeck mixed model. Both models are used in the comparative analysis of biological data originating from a set of living and/or extinct species to estimate the rate of phenotypic evolution resulting from genetic drift and selection. The function 'MiniBenchmark' allows to compare the calculation times for different tree sizes. See [this guide](https://venelin.github.io/SPLITT/articles/SPLITTRcppModules.html) for a tutorial on how to use the 'SPLITT' library in an R-package.

Installation
------------

You can install the released version of ThreePointUsingSPLITT from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ThreePointUsingSPLITT")
```

Examples
--------

-   Calculating the likelihood of the PMM for a given tree, data and model parameters.

``` r
library(ape)
library(ThreePointUsingSPLITT)
#> Loading required package: Rcpp

set.seed(10)

N <- 1000
x0 <- 0.1
alpha <- 1
theta <- 10
sigma2 <- 0.25
sigmae2 <- 1

tree <- rtree(N)

g <- rTraitCont(tree, model = "OU", root.value = x0,
                alpha = alpha, sigma = sqrt(sigma2),
                ancestor = FALSE)

x <- g + rnorm(n = N, mean = 0, sd = sqrt(sigmae2))

cat("logLikelihood using R:", POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2), "\n")
#> logLikelihood using R: -39915.58
cat("logLikelihood using R:", POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2), "\n")
#> logLikelihood using R: -39915.58
```

-   Performing a benchmark to measure the likelihood calculation times using different parallelization strategies:

``` r
# N specifies the size of the phylogenetic tree. 
# Ntests specifies the number of executions in one time measurement 
# (the more the better, but also slower).
MiniBenchmark(N = 100, Ntests = 100)
#> Performing a mini-benchmark of the POUMM log-likelihood calculation with 
#>       a tree of size N= 100 ;
#> Calling each likelihood calculation Ntests= 100  times ...
#> CPU:  (I don't know how to figure it out on your OS, sorry.) 
#> It seems that OpenMP was disabled during C++ compilation. For parallel
#>         tree traversal OpenMP should be available on your system and supported
#>         by your C++ compiler. Only serial tree traversal is possible with 
#>         the current SPLITT installation. Please, read the user guide for further 
#>         instructions on how to build the SPLITT library.Number of threads: 1 
#> Measuring calculation times...
#>    model                                            mode time.ms
#> 1  POUMM                                      R (serial)    1.10
#> 2  POUMM                                      C++ (AUTO)    0.03
#> 3  POUMM              C++ (SINGLE_THREAD_LOOP_POSTORDER)    0.02
#> 4  POUMM                 C++ (SINGLE_THREAD_LOOP_PRUNES)    0.03
#> 5  POUMM                 C++ (SINGLE_THREAD_LOOP_VISITS)    0.03
#> 6  POUMM                  C++ (MULTI_THREAD_LOOP_PRUNES)    0.03
#> 7  POUMM                  C++ (MULTI_THREAD_LOOP_VISITS)    0.02
#> 8  POUMM C++ (MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES)    0.03
#> 9  POUMM                  C++ (MULTI_THREAD_VISIT_QUEUE)    0.04
#> 10 POUMM     C++ (MULTI_THREAD_LOOP_PRUNES_NO_EXCEPTION)    0.02
#> 11 POUMM                        C++ (HYBRID_LOOP_PRUNES)    0.03
#> 12 POUMM                        C++ (HYBRID_LOOP_VISITS)    0.03
#> 13 POUMM       C++ (HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES)    0.03
MiniBenchmark(N = 1000, Ntests = 10)
#> Performing a mini-benchmark of the POUMM log-likelihood calculation with 
#>       a tree of size N= 1000 ;
#> Calling each likelihood calculation Ntests= 10  times ...
#> CPU:  (I don't know how to figure it out on your OS, sorry.) 
#> It seems that OpenMP was disabled during C++ compilation. For parallel
#>         tree traversal OpenMP should be available on your system and supported
#>         by your C++ compiler. Only serial tree traversal is possible with 
#>         the current SPLITT installation. Please, read the user guide for further 
#>         instructions on how to build the SPLITT library.Number of threads: 1 
#> Measuring calculation times...
#>    model                                            mode time.ms
#> 1  POUMM                                      R (serial)    12.0
#> 2  POUMM                                      C++ (AUTO)     0.2
#> 3  POUMM              C++ (SINGLE_THREAD_LOOP_POSTORDER)     0.2
#> 4  POUMM                 C++ (SINGLE_THREAD_LOOP_PRUNES)     0.1
#> 5  POUMM                 C++ (SINGLE_THREAD_LOOP_VISITS)     0.1
#> 6  POUMM                  C++ (MULTI_THREAD_LOOP_PRUNES)     0.1
#> 7  POUMM                  C++ (MULTI_THREAD_LOOP_VISITS)     0.1
#> 8  POUMM C++ (MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES)     0.1
#> 9  POUMM                  C++ (MULTI_THREAD_VISIT_QUEUE)     0.3
#> 10 POUMM     C++ (MULTI_THREAD_LOOP_PRUNES_NO_EXCEPTION)     0.1
#> 11 POUMM                        C++ (HYBRID_LOOP_PRUNES)     0.1
#> 12 POUMM                        C++ (HYBRID_LOOP_VISITS)     0.2
#> 13 POUMM       C++ (HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES)     0.1
MiniBenchmark(N = 10000, Ntests = 10)
#> Performing a mini-benchmark of the POUMM log-likelihood calculation with 
#>       a tree of size N= 10000 ;
#> Calling each likelihood calculation Ntests= 10  times ...
#> CPU:  (I don't know how to figure it out on your OS, sorry.) 
#> It seems that OpenMP was disabled during C++ compilation. For parallel
#>         tree traversal OpenMP should be available on your system and supported
#>         by your C++ compiler. Only serial tree traversal is possible with 
#>         the current SPLITT installation. Please, read the user guide for further 
#>         instructions on how to build the SPLITT library.Number of threads: 1 
#> Measuring calculation times...
#>    model                                            mode time.ms
#> 1  POUMM                                      R (serial)   126.0
#> 2  POUMM                                      C++ (AUTO)     1.4
#> 3  POUMM              C++ (SINGLE_THREAD_LOOP_POSTORDER)     1.6
#> 4  POUMM                 C++ (SINGLE_THREAD_LOOP_PRUNES)     1.5
#> 5  POUMM                 C++ (SINGLE_THREAD_LOOP_VISITS)     1.5
#> 6  POUMM                  C++ (MULTI_THREAD_LOOP_PRUNES)     1.6
#> 7  POUMM                  C++ (MULTI_THREAD_LOOP_VISITS)     1.9
#> 8  POUMM C++ (MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES)     1.4
#> 9  POUMM                  C++ (MULTI_THREAD_VISIT_QUEUE)     2.9
#> 10 POUMM     C++ (MULTI_THREAD_LOOP_PRUNES_NO_EXCEPTION)     1.3
#> 11 POUMM                        C++ (HYBRID_LOOP_PRUNES)     1.5
#> 12 POUMM                        C++ (HYBRID_LOOP_VISITS)     1.7
#> 13 POUMM       C++ (HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES)     1.5
```
