#include <Rcpp.h>
using namespace Rcpp;

// Integer power of 10 for exponents 0..9 (fits in 32-bit int)
static inline int pow10_int(int k) {
  int p = 1;
  for (int i = 0; i < k; ++i) p *= 10;
  return p;
}

// Count decimal digits of a non-negative int (0 has 1 digit)
static inline int ndigits_nonneg_int(int v) {
  if (v < 10) return 1;
  if (v < 100) return 2;
  if (v < 1000) return 3;
  if (v < 10000) return 4;
  if (v < 100000) return 5;
  if (v < 1000000) return 6;
  if (v < 10000000) return 7;
  if (v < 100000000) return 8;
  if (v < 1000000000) return 9;
  return 10; // up to INT_MAX has 10 digits
}

// [[Rcpp::export]]
List split_integer_first_last_matrix_cpp(IntegerMatrix x,
                                         int n_digits,
                                         int k_first,
                                         int k_last) {

  if (n_digits <= 0) stop("n_digits must be > 0.");
  if (k_first < 0 || k_last < 0) stop("k_first and k_last must be >= 0.");
  if (k_first + k_last > n_digits) stop("k_first + k_last must be <= n_digits.");

  const int exp_div = n_digits - k_first; // divisor exponent
  const int exp_mod = k_last;             // modulus exponent

  if (k_first > 0 && exp_div > 9)
    stop("10^(n_digits - k_first) would overflow 32-bit int; reduce n_digits or increase k_first.");
  if (k_last > 0 && exp_mod > 9)
    stop("10^(k_last) would overflow 32-bit int; reduce k_last.");

  const int div_first = (k_first == 0) ? 1 : pow10_int(exp_div);
  const int mod_last  = (k_last  == 0) ? 1 : pow10_int(exp_mod);

  const int nr = x.nrow();
  const int nc = x.ncol();

  IntegerMatrix first(nr, nc);
  IntegerMatrix last(nr, nc);

  // Column-major traversal (R matrix layout)
  for (int j = 0; j < nc; ++j) {
    for (int i = 0; i < nr; ++i) {
      int xi = x(i, j);

      if (IntegerMatrix::is_na(xi)) {
        first(i, j) = NA_INTEGER;
        last(i, j)  = NA_INTEGER;
        continue;
      }

      int ax = (xi < 0) ? -xi : xi;

      // Strict digit-length check
      int nd = ndigits_nonneg_int(ax);
      if (nd != n_digits) {
        stop("Element x[%d,%d] = %d has %d digits; expected exactly %d digits.",
             i + 1, j + 1, xi, nd, n_digits);
      }

      first(i, j) = (k_first == 0) ? 0 : ax / div_first;
      last(i, j)  = (k_last  == 0) ? 0 : ax % mod_last; // drops leading zeros in trailing group
    }
  }

  return List::create(_["first"] = first,
                      _["last"]  = last);
}
