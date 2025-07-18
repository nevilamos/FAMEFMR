// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List extract_season_firetypes(NumericMatrix v, CharacterVector colnames) {
  int nrows = v.nrow();
  int ncols = v.ncol();
  
  NumericMatrix seasons(nrows, ncols);
  NumericMatrix firetypes(nrows, ncols);
  
  for (int j = 0; j < ncols; ++j) {
    std::string name = Rcpp::as<std::string>(colnames[j]);
    
    int season = std::stoi(name.substr(0, 4));
    int firetype = std::stoi(name.substr(5, 1));  // 6th char, zero-indexed
    
    for (int i = 0; i < nrows; ++i) {
      double val = v(i, j);
      seasons(i, j) = val * season;
      firetypes(i, j) = val * firetype;
    }
  }
  
  return List::create(
    _["seasons"] = seasons,
    _["firetypes"] = firetypes
  );
}
