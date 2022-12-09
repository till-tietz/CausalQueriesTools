// [[Rcpp::depends(BH, bigmemory)]]
#include <RcppArmadillo.h>
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

// [[Rcpp::export]]
void realise_outcomes_c(SEXP outcomes,
                        std::vector<std::string> nodes,
                        std::vector<std::string> endogenous_nodes,
                        List dos,
                        List parents_list,
                        List nodal_types,
                        List nodal_types_colnames,
                        List nodal_types_collapsed,
                        int n_causal_types){

  //connect to bigmat
  Rcpp::XPtr<BigMatrix> outcomes_mat(outcomes);
  MatrixAccessor<int> outcomes_mat_access(*outcomes_mat);

  // get causal types
  std::vector<std::vector<std::string>> ct = make_causal_types_c(nodal_types_collapsed);
  // fill in values for dos
  std::vector<std::string> in_dos = dos.names();
  for(int i = 0; i < in_dos.size(); ++i){
    std::string in_dos_i = in_dos[i];
    int pos = std::find(nodes.begin(), nodes.end(), in_dos_i) - nodes.begin();
    int dos_i_int = dos[in_dos_i];
    std::vector<std::string> dos_i;
    dos_i.push_back(std::to_string(dos_i_int));
    ct[pos] = rep_times(dos_i, n_causal_types);
  }

  // loop over each endogenous node
  for(int i = 0; i < endogenous_nodes.size(); ++i){
    std::string var = endogenous_nodes[i];
    int pos = std::find(nodes.begin(), nodes.end(), var) - nodes.begin();
    //get causal type realizations for endogenous node
    std::vector<std::string> child_type = ct[pos];
    //get parents of endogenous node
    std::vector<std::string> parents = parents_list[var];
    //get nodal types of endogenous node
    std::vector<std::string> nodal_label = nodal_types_collapsed[var];
    //get uncollapsed nodal types for endogenous node
    arma::mat nodal_type_var = nodal_types[var];
    //get uncollapsed nodal types colnames
    std::vector<std::string> nodal_type_var_col = nodal_types_colnames[var];

    //loop over causal types
    for(int j = 0; j < child_type.size(); ++j){
      //get causal type
      std::string type = child_type[j];
      //generate empty vector for parent realization
      std::string parents_val;
      //get parent realization
      for(int k = 0; k < parents.size(); ++k){
        int pos_parent = std::find(nodes.begin(), nodes.end(), parents[k]) - nodes.begin();
        parents_val.append(ct[pos_parent][j]);
      }
      //find row position of type
      int row = std::find(nodal_label.begin(),nodal_label.end(), type) - nodal_label.begin();
      //find realization and add to J
      int pos_col = std::find(nodal_type_var_col.begin(), nodal_type_var_col.end(), parents_val) - nodal_type_var_col.begin();
      int outcome = int(nodal_type_var(row,pos_col));
      ct[pos][j] = std::to_string(outcome);
    }
  }

  for(int i = 0; i < ct.size(); ++i){
    for(int j = 0; j < n_causal_types; ++j){
      outcomes_mat_access[i][j] = std::stoi(ct[i][j]);
    }
  }

  return;
}


//[[Rcpp::export]]
void realise_outcomes_singular_c(SEXP outcomes,
                                 std::vector<std::string> nodes,
                                 std::vector<std::string> endogenous_nodes,
                                 List dos,
                                 List parents_list,
                                 List nodal_types,
                                 List nodal_types_colnames,
                                 List nodal_types_collapsed,
                                 int n_causal_types,
                                 std::vector<std::string> vars){


 //connect to bigmat
 Rcpp::XPtr<BigMatrix> outcomes_mat(outcomes);
 MatrixAccessor<int> outcomes_mat_access(*outcomes_mat);

 std::vector<std::string> in_dos = dos.names();

 for(int i = 0; i < dos.size(); ++i){
   // generate causal types and fill in dos
   std::vector<std::vector<std::string>> ct = make_causal_types_c(nodal_types_collapsed);
   std::string in_dos_i = in_dos[i];
   int pos = std::find(nodes.begin(), nodes.end(), in_dos_i) - nodes.begin();
   int dos_i_int = dos[in_dos_i];
   std::vector<std::string> dos_i;
   dos_i.push_back(std::to_string(dos_i_int));
   ct[pos] = rep_times(dos_i, n_causal_types);

   // loop over each endogenous node
   for(int j = 0; j < endogenous_nodes.size(); ++j){
     std::string var = endogenous_nodes[j];
     int pos = std::find(nodes.begin(), nodes.end(), var) - nodes.begin();
     //get causal type realizations for endogenous node
     std::vector<std::string> child_type = ct[pos];
     //get parents of endogenous node
     std::vector<std::string> parents = parents_list[var];
     //get nodal types of endogenous node
     std::vector<std::string> nodal_label = nodal_types_collapsed[var];
     //get uncollapsed nodal types for endogenous node
     arma::mat nodal_type_var = nodal_types[var];
     //get uncollapsed nodal types colnames
     std::vector<std::string> nodal_type_var_col = nodal_types_colnames[var];

     //loop over causal types
     for(int k = 0; k < child_type.size(); ++k){
       //get causal type
       std::string type = child_type[k];
       //generate empty vector for parent realization
       std::string parents_val;
       //get parent realization
       for(int l = 0; l < parents.size(); ++l){
         int pos_parent = std::find(nodes.begin(), nodes.end(), parents[l]) - nodes.begin();
         parents_val.append(ct[pos_parent][k]);
       }
       //find row position of type
       int row = std::find(nodal_label.begin(),nodal_label.end(), type) - nodal_label.begin();
       //find realization and add to J
       int pos_col = std::find(nodal_type_var_col.begin(), nodal_type_var_col.end(), parents_val) - nodal_type_var_col.begin();
       int outcome = int(nodal_type_var(row,pos_col));
       ct[pos][k] = std::to_string(outcome);
     }
   }

   std::string var_i = vars[i];
   pos = std::find(nodes.begin(), nodes.end(), var_i) - nodes.begin();

   for(int m = 0; m < n_causal_types; ++m){
     outcomes_mat_access[i][m] = std::stoi(ct[pos][m]);
   }
 }

 return;

}

