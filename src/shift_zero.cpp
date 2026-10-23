#include <Rcpp.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void shift_zero_in_place_cpp(IntegerMatrix x)
{
  int m = x.nrow();
  int n = x.ncol();

  for (int i = 0, k = 0, k0 = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {

      // Keep only non-zero values (zeros treated as empty)
      if (x[k] != 0) {
        x[k0] = x[k];
        k0 += m;
      }
      k += m;
    }

    // Fill remainder with 0
    while (k0 < k) {
      x[k0] = 0;
      k0 += m;
    }

    k  = (k % m) + 1;
    k0 = k;
  }

  // Drop column names (dimnames[[2]]) as in original
  if (x.attr("dimnames") != R_NilValue) {
    List dn = x.attr("dimnames");
    dn[1] = R_NilValue;

    if (dn.attr("names") != R_NilValue) {
      CharacterVector ndn = dn.attr("names");
      ndn[1] = "";
    }
  }
}

// [[Rcpp::export]]
IntegerMatrix shift_zero_cpp(IntegerMatrix x)
{
  IntegerMatrix y = clone(x);
  shift_zero_in_place_cpp(y);
  return y;
}
