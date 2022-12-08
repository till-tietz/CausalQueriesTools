// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>
using namespace Rcpp;

//' cpp implementation of standard R rep
//'
//' @param x integer vector to be repeated
//' @param n integer number of repetitions
//' @return integer vector x repeated n times
//' @keywords internal
// [[Rcpp::export]]
std::vector<std::string> rep_times(std::vector<std::string> x, int n){

  std::vector<std::string> ret;

  for(int i = 0; i < n; ++i){
    ret.insert(ret.end(),x.begin(),x.end());
  }

  return ret;
}

//' cpp implementation of R rep with each argument
//'
//' @param x char vector of elements to be repeated
//' @param n integer number of repetitions
//' @return char vector with each element in x repeated n times in order
//' @keywords internal
// [[Rcpp::export]]
std::vector<std::string> rep_each(std::vector<std::string> x, int n){

  std::vector<std::string> ret(x.size() * n);
  int ind = -1;

  for(int i = 0; i < x.size(); ++i){
    for(int j = 0; j < n; ++j){
      ind += 1;
      ret[ind] = x[i];
    }
  }

  return ret;
}

//' cpp helper to make causal types
//'
//' @param nodal_types a List of nodal types
//' @return vector of vectors containing causal types
//' @keywords internal
// [[Rcpp::export]]
std::vector<std::vector<std::string>> make_causal_types_c(List nodal_types){

  //number of nodal types per node
  std::vector<int> n_nodal_types(nodal_types.size());

  for(int i = 0; i < nodal_types.size(); ++i){
    std::vector<std::string> nodal_types_i = nodal_types[i];
    n_nodal_types[i] = nodal_types_i.size();
  }

  n_nodal_types.insert(n_nodal_types.begin(),1);
  n_nodal_types.insert(n_nodal_types.end(),1);

  std::vector<std::vector<std::string>> causal_types(nodal_types.size());

  for(int i = n_nodal_types.size() - 2; i > 0; --i){
    std::vector<std::string> ct_i;
    int each = std::accumulate(std::begin(n_nodal_types), std::begin(n_nodal_types) + i, 1, std::multiplies<int>());
    ct_i = rep_each(nodal_types[i-1],each);
    int times = std::accumulate(std::begin(n_nodal_types) + i + 1, std::end(n_nodal_types), 1, std::multiplies<int>());
    ct_i = rep_times(ct_i,times);
    causal_types[i-1] = ct_i;
  }

  return causal_types;
}

