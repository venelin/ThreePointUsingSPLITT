#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]

// BEGIN: Needed for r-devel (R 3.4)
void R_init_ThreePointUsingSPLITT(DllInfo *info) {
  /* Register routines, allocate resources. */
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}


void R_unload_ThreePointUsingSPLITT(DllInfo *info) {
  /* Release resources. */
}
// END Needed for r-devel (R 3.4)
