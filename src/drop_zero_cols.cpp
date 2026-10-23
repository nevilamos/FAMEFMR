#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix drop_zero_cols_cpp(const IntegerMatrix& mat) {
  const int nRows = mat.nrow();
  const int nCols = mat.ncol();

  std::vector<bool> keepCol(nCols, false);
  int nKeep = 0;

  // 1) Identify columns that are not all zero
  for (int j = 0; j < nCols; ++j) {
    for (int i = 0; i < nRows; ++i) {
      if (mat(i, j) != 0) {
        keepCol[j] = true;
        ++nKeep;
        break;
      }
    }
  }

  // 2) Create new matrix and fill it
  IntegerMatrix res(nRows, nKeep);
  int currentCol = 0;

  for (int j = 0; j < nCols; ++j) {
    if (keepCol[j]) {
      res(_, currentCol) = mat(_, j);
      ++currentCol;
    }
  }

  // Optional: preserve rownames; drop colnames (or preserve selected colnames)
  if (mat.hasAttribute("dimnames")) {
    List dn = mat.attr("dimnames");
    // Keep rownames if present
    if (dn.size() >= 1) {
      List dn2(2);
      dn2[0] = dn[0];        // rownames
      dn2[1] = R_NilValue;   // colnames not preserved by default
      res.attr("dimnames") = dn2;
    }
  }

  return res;
}
