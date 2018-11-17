library(testthat)
context("Test R and Cpp code calculate the same PMM log-likelihood")

library(ape)
library(PMMUsingSPLITT)

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

cppObj3Point <- New3PointPOUMMCppObject(x, tree)

test_that(
  "POUMMLogLik == POUMMLogLikCpp 3-point", {
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObj3Point, 0))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObj3Point, 10))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObj3Point, 11))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObj3Point, 12))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObj3Point, 21))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObj3Point, 22))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObj3Point, 23))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObj3Point, 24))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObj3Point, 25))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObj3Point, 31))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObj3Point, 32))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObj3Point, 33))
  })

cppObjAbc <- NewAbcPOUMMCppObject(x, tree)
test_that(
  "POUMMLogLik == POUMMLogLikCpp Abc", {
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObjAbc, 0))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObjAbc, 10))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObjAbc, 11))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObjAbc, 12))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObjAbc, 21))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObjAbc, 22))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObjAbc, 23))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObjAbc, 24))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObjAbc, 25))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObjAbc, 31))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObjAbc, 32))
    expect_equal(POUMMLogLik(x, tree, x0, alpha, theta, sigma2, sigmae2),
                 POUMMLogLikCpp(x, tree, x0, alpha, theta, sigma2, sigmae2, 
                                cppObjAbc, 33))
  }
)
