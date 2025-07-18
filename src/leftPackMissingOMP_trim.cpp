// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix leftPackMissingOMP_trim(NumericMatrix mat) {
  int nrows = mat.nrow();
  int ncols = mat.ncol();
  
  // Step 1: Count non-NA per row & find max
  IntegerVector row_counts(nrows);
  int max_cols = 0;
  
#pragma omp parallel for reduction(max:max_cols)
  for (int i = 0; i < nrows; ++i) {
    int count = 0;
    for (int j = 0; j < ncols; ++j) {
      if (!NumericVector::is_na(mat(i, j))) {
        count++;
      }
    }
    row_counts[i] = count;
    if (count > max_cols) max_cols = count;
  }
  
  // Step 2: Create output matrix filled with NA
  NumericMatrix out(nrows, max_cols);
  std::fill(out.begin(), out.end(), NA_REAL);
  
  // Step 3: Left-pack non-NA values
#pragma omp parallel for
  for (int i = 0; i < nrows; ++i) {
    int pos = 0;
    for (int j = 0; j < ncols; ++j) {
      double val = mat(i, j);
      if (!NumericVector::is_na(val)) {
        out(i, pos++) = val;
      }
    }
  }
  
  return out;
}
