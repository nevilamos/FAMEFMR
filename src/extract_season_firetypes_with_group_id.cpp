// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

// Hash a row to a string key
std::string row_key(const NumericVector& row) {
  std::string key;
  for (int i = 0; i < row.size(); ++i) {
    if (NumericVector::is_na(row[i]))
      key += "NA,";
    else
      key += std::to_string(row[i]) + ",";
  }
  return key;
}

// [[Rcpp::export]]
List extract_season_firetypes_with_group_id(NumericMatrix v, CharacterVector colnames) {
  int nrows = v.nrow();
  int ncols = v.ncol();
  
  NumericMatrix seasons(nrows, ncols);
  NumericMatrix firetypes(nrows, ncols);
  IntegerVector group_ids(nrows);
  
  // Step 1: Extract season and firetype matrices
  for (int j = 0; j < ncols; ++j) {
    std::string name = Rcpp::as<std::string>(colnames[j]);
    int season = std::stoi(name.substr(0, 4));
    int firetype = std::stoi(name.substr(5, 1));
    
    for (int i = 0; i < nrows; ++i) {
      double val = v(i, j);
      seasons(i, j) = val * season;
      firetypes(i, j) = val * firetype;
    }
  }
  
  // Step 2: Assign group_id by hashing the season matrix rows
  std::unordered_map<std::string, int> key_to_group;
  int group_counter = 1;
  
  for (int i = 0; i < nrows; ++i) {
    NumericVector row = seasons(i, _);
    std::string key = row_key(row);
    
    auto it = key_to_group.find(key);
    if (it == key_to_group.end()) {
      key_to_group[key] = group_counter;
      group_ids[i] = group_counter;
      group_counter++;
    } else {
      group_ids[i] = it->second;
    }
  }
  
  return List::create(
    _["FireTypeNos"] = firetypes,
    _["SEASONS"] = seasons,
    _["group_id"] = group_ids
  );
}
