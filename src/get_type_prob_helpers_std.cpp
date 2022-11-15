#include <Rcpp.h>

using namespace Rcpp;

//' generates one draw from type probability distribution for each type in P
//'
//' @param P parameter_matrix of parameters and causal types
//' @param parameters, priors or posteriors
//' @param nrow number of rows in P
//' @param ncol number of columns in P
//' @return draw from type distribution for each type in P
// [[Rcpp::export]]
std::vector<double> get_type_prob_cstd(std::vector<double> P, std::vector<double> parameters, int nrow, int ncol){

  for (int j=0; j<ncol; j++) {
    for (int i=0; i<nrow; i++)  {
      int id = (j+1) * nrow - (nrow - i);
      P[id] = P[id] * parameters[i] + 1 - P[id];
    }
  }

  std::vector<double> ret;

  for (int j=0; j<ncol; j++) {
    int start = (j+1) * nrow - nrow;
    int end = (j+1) * nrow;
    std::vector<double> v_j(P.begin() + start, P.begin() + end);
    ret.push_back(std::accumulate(v_j.begin(), v_j.end(), 1.0, std::multiplies<double>()));
  }
  return ret;
}

//' generates n draws from type probability distribution for each type in P
//'
//' @param params parameters, priors or posteriors
//' @param P parameter_matrix of parameters and causal types
//' @param nrow number of rows in params matrix before being flattend to vector
//' @param ncol numbers of columns in params matrix before being flattened to vector
//' @param nrow_p number of rows in P matrix before being flattened to vector
//' @param ncol_p number of columns in P matrix befroe being flattened to vector
//' @return draws from type distribution for each type in P
// [[Rcpp::export]]
std::vector<double> get_type_prob_multiple_cstd(std::vector<double> params, std::vector<double> P, int nrow, int ncol, int nrow_p, int ncol_p){

  std::vector<double> ret;

  for (int i=0; i<nrow; i++){

    std::vector<double> row_i;

    for (int j = 0; j < ncol; j++){
      row_i.push_back(params[j * nrow + i]);
    }

    std::vector<double> col_i = get_type_prob_cstd(P,row_i,ncol_p,nrow_p);
    std::copy(col_i.begin(), col_i.end(), std::back_inserter(ret));
  }

  return ret;
}
