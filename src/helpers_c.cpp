// [[Rcpp::depends(BH, bigmemory)]]
#include <RcppArmadillo.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>
using namespace Rcpp;


//***************************
//* generating causal types *
//***************************

// cpp implementation of standard R rep
std::vector<std::string> rep_times(std::vector<std::string> x, int n){

  std::vector<std::string> ret;

  for(int i = 0; i < n; ++i){
    ret.insert(ret.end(),x.begin(),x.end());
  }

  return ret;
}

// cpp implementation of R rep with each argument
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

// cpp helper to make causal types
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

//********************
//* realise outcomes *
//********************

//' cpp implementation of realise_outcomes. Realise outcomes for all causal types.
//' Calculated by sequentially calculating endogenous nodes. If a do operator is applied to
//' any node then it takes the given value and all its descendants are generated accordingly.
//' Output is written to a bigmatrix.
//'
//' @param outcomes memory address of a bigmatrix object
//' @param nodes string vector of nodes names
//' @param endogenous_nodes string vector of endogenous nodes
//' @param dos List of do operations
//' @param parents_list List of parent nodes for each node
//' @param nodal_types List of integer matrices with uncollapsed nodal types
//' @param nodal_types_colnames List of column names of matrices in nodal_types
//' @param nodal_types_collapsed List of collapsed nodal types
//' @param n_causal_types int specifying number of causal types
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

//*******************************
//* map queries to causal types *
//*******************************

// helper: element wise vector addition
std::vector<int> add(std::vector<int> x,
                     std::vector<int> y){

  if(y.size() == 1){
    for(int i = 0; i < x.size(); ++i){
      x[i] += y[0];
    }
  } else {
    for(int i = 0; i < x.size(); ++i){
      x[i] += y[i];
    }
  }
  return x;
}

// helper: element wise vector subtraction
std::vector<int> subtract(std::vector<int> x,
                          std::vector<int> y){

  if(y.size() == 1){
    for(int i = 0; i < x.size(); ++i){
      x[i] -= y[0];
    }
  } else {
    for(int i = 0; i < x.size(); ++i){
      x[i] -= y[i];
    }
  }
  return x;
}

// helper: element wise strictly greater than
std::vector<int> strict_greater_than(std::vector<int> x,
                                     std::vector<int> y){

  if(y.size() == 1){
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] > y[0];
      x[i] = x_i;
    }
  } else {
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] > y[i];
      x[i] = x_i;
    }
  }
  return x;
}

// helper: element wise greater than
std::vector<int> greater_than(std::vector<int> x,
                              std::vector<int> y){
  if(y.size() == 1){
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] >= y[0];
      x[i] = x_i;
    }
  } else {
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] >= y[i];
      x[i] = x_i;
    }
  }
  return x;
}

// helper: element wise strictly smaller than
std::vector<int> strict_smaller_than(std::vector<int> x,
                                     std::vector<int> y){
  if(y.size() == 1){
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] < y[0];
      x[i] = x_i;
    }
  } else {
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] < y[i];
      x[i] = x_i;
    }
  }
  return x;
}

// helper: element wise smaller than
std::vector<int> smaller_than(std::vector<int> x,
                              std::vector<int> y){
  if(y.size() == 1){
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] <= y[0];
      x[i] = x_i;
    }
  } else {
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] <= y[i];
      x[i] = x_i;
    }
  }
  return x;
}

// helper: element wise equal
std::vector<int> equal(std::vector<int> x,
                       std::vector<int> y){
  if(y.size() == 1){
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] == y[0];
      x[i] = x_i;
    }
  } else {
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] == y[i];
      x[i] = x_i;
    }
  }
  return x;
}

// helper: element wise not equal
std::vector<int> not_equal(std::vector<int> x,
                           std::vector<int> y){
  if(y.size() == 1){
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] != y[0];
      x[i] = x_i;
    }
  } else {
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] != y[i];
      x[i] = x_i;
    }
  }
  return x;
}

// helper: element wise and
std::vector<int> And(std::vector<int> x,
                     std::vector<int> y){
  if(y.size() == 1){
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] & y[0];
      x[i] = x_i;
    }
  } else {
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] & y[i];
      x[i] = x_i;
    }
  }
  return x;
}

//helper: element wise or
std::vector<int> Or(std::vector<int> x,
                    std::vector<int> y){
  if(y.size() == 1){
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] | y[0];
      x[i] = x_i;
    }
  } else {
    for(int i = 0; i < x.size(); ++i){
      int x_i = x[i] | y[i];
      x[i] = x_i;
    }
  }
  return x;
}


std::vector<int> pair_operation(std::vector<int> a,
                                std::vector<int> b,
                                std::string op){
  std::vector<int> ret;

  if(op == "+"){
    ret = add(a,b);
  }

  if(op == "-"){
    ret = subtract(a,b);
  }

  if(op == ">"){
    ret = strict_greater_than(a,b);
  }

  if(op == ">="){
    ret = greater_than(a,b);
  }

  if(op == "<"){
    ret = strict_smaller_than(a,b);
  }

  if(op == "<="){
    ret = smaller_than(a,b);
  }

  if(op == "=="){
    ret = equal(a,b);
  }

  if(op == "!="){
    ret = not_equal(a,b);
  }

  if(op == "&" || op == "&&"){
    ret = And(a,b);
  }

  if(op == "|" || op == "||"){
    ret = Or(a,b);
  }
  return ret;
}


//' cpp implementation of realise_outcomes for map_query_to_causal_types. Dos are evaluated
//' and the realised outcomes for the variable they are attached to is written to a bigmatrix.
//'
//' @param nodes string vector of nodes names
//' @param endogenous_nodes string vector of endogenous nodes
//' @param dos List of do operations
//' @param parents_list List of parent nodes for each node
//' @param nodal_types List of integer matrices with uncollapsed nodal types
//' @param nodal_types_colnames List of column names of matrices in nodal_types
//' @param nodal_types_collapsed List of collapsed nodal types
//' @param n_causal_types int specifying number of causal types
//' @param vars string vector with names of variables dos are attached to
// [[Rcpp::export]]
std::vector<std::vector<int>> query_to_ct_c(
                                 std::vector<std::string> nodes,
                                 std::vector<std::string> endogenous_nodes,
                                 List dos,
                                 List parents_list,
                                 List nodal_types,
                                 List nodal_types_colnames,
                                 List nodal_types_collapsed,
                                 int n_causal_types,
                                 std::vector<std::string> vars){

  //vector of vectors to store data realisations
  std::vector<std::vector<int>> out_mat;
  for(int i = 0; i < vars.size(); ++i){
    out_mat.push_back(std::vector<int> (n_causal_types));
  }

  std::vector<std::string> in_dos = dos.names();

  for(int i = 0; i < dos.size(); ++i){
    // generate causal types and fill in dos
    std::vector<std::vector<std::string>> ct = make_causal_types_c(nodal_types_collapsed);
    std::string in_dos_i = in_dos[i];
    int pos = std::find(nodes.begin(), nodes.end(), in_dos_i) - nodes.begin();
    int dos_i_int = dos[i];

    if(dos_i_int >= 0){
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
    }

    std::string var_i = vars[i];
    pos = std::find(nodes.begin(), nodes.end(), var_i) - nodes.begin();
    //put realzations into out_mat
    for(int m = 0; m < n_causal_types; ++m){
      out_mat[i][m] = std::stoi(ct[pos][m]);
    }
  }

  return out_mat;

}

