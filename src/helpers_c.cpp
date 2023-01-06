// [[Rcpp::depends(BH, bigmemory)]]
#include <RcppArmadillo.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>
using namespace Rcpp;


//***************************
//* generating causal types *
//***************************

// cpp implementation of standard R rep
//std::vector<std::string> rep_times(std::vector<std::string> x, int n){
//
//  std::vector<std::string> ret;

//  for(int i = 0; i < n; ++i){
//    ret.insert(ret.end(),x.begin(),x.end());
//  }
//
//  return ret;
//}

// cpp implementation of R rep with each argument
//std::vector<std::string> rep_each(std::vector<std::string> x, int n){

//  std::vector<std::string> ret(x.size() * n);
//  int ind = -1;
//
//  for(int i = 0; i < x.size(); ++i){
//    for(int j = 0; j < n; ++j){
//      ind += 1;
//      ret[ind] = x[i];
//    }
//  }
//
//  return ret;
//}

// cpp helper to make causal types
//void make_causal_types_c_test(SEXP ct_mat_address,
//                              List nodal_types,
//                              std::vector<std::string> unique_nodal_types){
//
  //connect to bigmat
//  Rcpp::XPtr<BigMatrix> ct_mat(ct_mat_address);
//  MatrixAccessor<int> ct_mat_access(*ct_mat);

  //number of nodal types per node
//  std::vector<int> n_nodal_types(nodal_types.size());
//
//  for(int i = 0; i < nodal_types.size(); ++i){
//    std::vector<std::string> nodal_types_i = nodal_types[i];
//    n_nodal_types[i] = nodal_types_i.size();
//  }
//
//  n_nodal_types.insert(n_nodal_types.begin(),1);
//  n_nodal_types.insert(n_nodal_types.end(),1);
//
//  for(int i = n_nodal_types.size() - 2; i > 0; --i){
//
//    int each = std::accumulate(std::begin(n_nodal_types), std::begin(n_nodal_types) + i, 1, std::multiplies<int>());
//    int times = std::accumulate(std::begin(n_nodal_types) + i + 1, std::end(n_nodal_types), 1, std::multiplies<int>());
//    std::vector<std::string> nodal_types_i = nodal_types[i-1];
//    std::vector<int> nodal_types_int(nodal_types_i.size());
//
//    for(int j = 0; j < nodal_types_i.size(); ++j){
//      nodal_types_int[j] = std::find(unique_nodal_types.begin(), unique_nodal_types.end(), nodal_types_i[j]) - unique_nodal_types.begin();
//    }
//
//    int pos = 0;
//    for(int j = 0; j < times; ++j){
//      for(int k = 0; k < nodal_types_int.size(); ++k){
//        for(int l = 0; l < each; ++l){
//          ct_mat_access[i-1][pos] = nodal_types_int[k];
//          pos += 1;
//        }
//      }
//    }
//
//  }
//  return;
//}


//' Determines the number of times each nodes' vector of nodal types is repeated to generate the pattern of causal types.
//' Combining this information with the position of a causal type allows us to determine the nodal types used in the construction
//' of the given causal type.
//'
//' @param nodal_types a List of nodal types
//' @return an integer vector of repetitions for each nodes' nodal type vector
// [[Rcpp::export]]
std::vector<int> get_causal_type_pattern(List nodal_types){

  //number of nodal types per node
  std::vector<int> n_nodal_types(nodal_types.size());

  for(int i = 0; i < nodal_types.size(); ++i){
    std::vector<std::string> nodal_types_i = nodal_types[i];
    n_nodal_types[i] = nodal_types_i.size();
  }

  n_nodal_types.insert(n_nodal_types.begin(),1);
  n_nodal_types.insert(n_nodal_types.end(),1);
  std::vector<int> ret;

  for(int i = n_nodal_types.size() - 2; i > 0; --i){
    ret.push_back(std::accumulate(std::begin(n_nodal_types), std::begin(n_nodal_types) + i, 1, std::multiplies<int>()));
  }

  return ret;
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
                        SEXP causal_types,
                        std::vector<std::string> nodes,
                        std::vector<std::string> endogenous_nodes,
                        List dos,
                        List parents_list,
                        List nodal_types,
                        List nodal_types_colnames,
                        List nodal_types_collapsed,
                        std::vector<std::string> unique_nodal_types,
                        int n_causal_types){

  //connect to output bigmat
  Rcpp::XPtr<BigMatrix> outcomes_mat(outcomes);
  MatrixAccessor<int> outcomes_mat_access(*outcomes_mat);

  //connect to causal types bigmat
  Rcpp::XPtr<BigMatrix> ct_mat(causal_types);
  MatrixAccessor<int> ct(*ct_mat);

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
//' @param List of pairwise operations query is made up of
// [[Rcpp::export]]
void query_to_ct_c(SEXP outcomes,
                   std::vector<std::string> nodes,
                   std::vector<std::string> endogenous_nodes,
                   List dos,
                   List parents_list,
                   List nodal_types,
                   List nodal_types_colnames,
                   List nodal_types_collapsed,
                   std::vector<std::string> unique_nodal_types,
                   int n_causal_types,
                   std::vector<std::string> vars,
                   List query_operations,
                   std::vector<std::string> var_order){
  //connect to bigmat
  Rcpp::XPtr<BigMatrix> outcomes_mat(outcomes);
  MatrixAccessor<int> outcomes_mat_access(*outcomes_mat);

  for(int i = 0; i < dos.size(); ++i){
    // generate causal types and fill in dos
    List dos_i = dos[i];
    std::vector<std::string> in_dos_i = dos_i.names();

    if(dos_i.size() == 1){
      int dos_i_val = dos_i[0];
      if(dos_i_val < 0){
        in_dos_i[0] = "";
      }
    }

    std::vector<std::string> work_through;
    for(int n = 0; n < endogenous_nodes.size(); ++n){
      for(int o = 0; o < in_dos_i.size(); ++o){
        if(endogenous_nodes[n] != in_dos_i[o]){
          work_through.push_back(endogenous_nodes[n]);
        }
      }
    }


    for(int j = 0; j < n_causal_types; ++j){

      for(int k = 0; k < work_through.size(); ++k){
        std::string var = work_through[k];
        int pos = std::find(nodes.begin(), nodes.end(), var) - nodes.begin();
      }

    }


    // loop over each endogenous node
    for(int j = 0; j < work_through.size(); ++j){
      std::string var = work_through[j];
      int pos = std::find(nodes.begin(), nodes.end(), var) - nodes.begin();
      //get parents of endogenous node
      std::vector<std::string> parents = parents_list[var];
      //get nodal types of endogenous node
      std::vector<std::string> nodal_label = nodal_types_collapsed[var];
      //get uncollapsed nodal types for endogenous node
      arma::mat nodal_type_var = nodal_types[var];
      //get uncollapsed nodal types colnames
      std::vector<std::string> nodal_type_var_col = nodal_types_colnames[var];

      //loop over causal types
      for(int k = 0; k < n_causal_types; ++k){
        //get causal type
        int type_int = ct[pos][k];
        std::string type = unique_nodal_types[type_int];
        //generate empty vector for parent realization
        std::string parents_val;
        //get parent realization
        for(int l = 0; l < parents.size(); ++l){
          int parent_do = std::find(in_dos_i.begin(),in_dos_i.end(),parents[l]) - in_dos_i.begin();

          if(parent_do != in_dos_i.size()){
            std::string =
          }
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
    //put realzations into out_mat
    for(int m = 0; m < n_causal_types; ++m){
      out_mat[i][m] = std::stoi(ct[pos][m]);
    }
  }

  //apply mathematical operations
  //initial operations
  std::vector<std::vector<int>> op_init;

  //do initial individual binary arithmetic operations
  for(int i = 0; i < query_operations.size(); ++i){
    std::vector<std::string> op_i = query_operations[i];

    int pos_elem_1 = std::find(var_order.begin(), var_order.end(), op_i[0]) - var_order.begin();
    std::vector<int> elem_1;
    if(pos_elem_1 == var_order.size()){
      elem_1.push_back(std::stoi(op_i[0]));
    } else {
      elem_1 = out_mat[pos_elem_1];
    }

    int pos_elem_2 = std::find(var_order.begin(), var_order.end(), op_i[2]) - var_order.begin();
    std::vector<int> elem_2;
    if(pos_elem_2 == var_order.size()){
      elem_2.push_back(std::stoi(op_i[2]));
    } else {
      elem_2 = out_mat[pos_elem_2];
    }

    op_init.push_back(pair_operation(elem_1, elem_2, op_i[1]));

  }

  //combine binary arithmetic results with &
  std::vector<int> ret;

  if(op_init.size() > 1){
    //to do >> write & operations over op_init
    ret = pair_operation(op_init[0],op_init[1], "&");
    for(int i = 2; i < op_init.size(); ++i){
      ret = pair_operation(ret,op_init[i], "&");
    }
  } else {
    //if only one expression supplied simply return result of that operation
    ret = op_init[0];
  }

  //write to bigmatrix
  for(int i = 0; i < ret.size(); ++i){
    outcomes_mat_access[0][i] = ret[i];
  }

  return;

}

