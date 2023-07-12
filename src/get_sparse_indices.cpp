#include <Rcpp.h>
#include <vector>

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::List get_nonzero_row_indices_cpp(const std::vector<int>& idx_i, const std::vector<int>& idx_j) {
  int n = *std::max_element(idx_j.begin(), idx_j.end());
  std::vector<std::vector<int>> nz_indices(n);
  
  for (size_t k = 0; k < idx_i.size(); ++k) {
    int col = idx_j[k] - 1; // Convert 1-based to 0-based indexing
    int row = idx_i[k] - 1; // Convert 1-based to 0-based indexing
    nz_indices[col].push_back(row);
  }
  
  Rcpp::List result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = nz_indices[i];
  }
  
  return result;
}

// [[Rcpp::export]]
Rcpp::List get_nonzero_row_values_cpp(const std::vector<int>& y, const std::vector<int>& idx_j) {
  int n = *std::max_element(idx_j.begin(), idx_j.end());
  std::vector<std::vector<int>> nz_values(n);
  
  for (size_t k = 0; k < y.size(); ++k) {
    int col = idx_j[k] - 1; // Convert 1-based to 0-based indexing
    int y_val = y[k]; // Convert 1-based to 0-based indexing
    nz_values[col].push_back(y_val);
  }
  
  Rcpp::List result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = nz_values[i];
  }
  
  return result;
}
