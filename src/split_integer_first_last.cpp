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
  return 10; // up to INT_MAX (2,147,483,647)
}

// [[Rcpp::export]]
Rcpp::List split_integer_first_last_cpp(Rcpp::IntegerVector x,
                                        int n_digits,
                                        int k_first,
                                        int k_last) {

  if (n_digits <= 0) stop("n_digits must be > 0.");
  if (k_first < 0 || k_last < 0) stop("k_first and k_last must be >= 0.");
  if (k_first + k_last > n_digits)
    stop("k_first + k_last must be <= n_digits.");

  // powers of 10 must fit in 32-bit int
  int exp_div = n_digits - k_first;
  int exp_mod = k_last;

  if (k_first > 0 && exp_div > 9)
    stop("10^(n_digits - k_first) would overflow 32-bit int.");
  if (k_last > 0 && exp_mod > 9)
    stop("10^(k_last) would overflow 32-bit int.");

  const int div_first = (k_first == 0) ? 1 : pow10_int(exp_div);
  const int mod_last  = (k_last  == 0) ? 1 : pow10_int(exp_mod);

  int N = x.size();
  IntegerVector first(N);
  IntegerVector last(N);

  for (int i = 0; i < N; ++i) {
    int xi = x[i];

    if (IntegerVector::is_na(xi)) {
      first[i] = NA_INTEGER;
      last[i]  = NA_INTEGER;
      continue;
    }

    // absolute value for digit counting and splitting
    int ax = (xi < 0) ? -xi : xi;

    // strict digit-count check
    int nd = ndigits_nonneg_int(ax);
    if (nd != n_digits) {
      stop("Element x[%d] = %d has %d digits; expected exactly %d digits.",
           i + 1, xi, nd, n_digits);
    }

    first[i] = (k_first == 0) ? 0 : ax / div_first;
    last[i]  = (k_last  == 0) ? 0 : ax % mod_last;
  }

  return List::create(
    _["first"] = first,
    _["last"]  = last
  );
}
