// [[Rcpp::depends(BH, bigmemory)]]
#include <RcppArmadillo.h>
#include <numeric>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
void query_to_ct_c(std::vector<std::string> nodes,
                   List nodal_types,
                   int n_causal_types){

  int n_nodes = nodes.size();
  for(int i = 0; i < n_causal_types; ++i){

    //generate causal type components ---------------------------

    //number of nodal types per node
    std::vector<int> n_nodal_types(n_nodes);

    for(int j = 0; j < n_nodes; ++j){
      std::vector<std::string> nodal_types_j = nodal_types[j];
      n_nodal_types[j] = nodal_types_j.size();
    }

    //number of times nodal types are repeated
    n_nodal_types.insert(n_nodal_types.begin(),1);
    n_nodal_types.insert(n_nodal_types.end(),1);
    std::vector<int> reps;

    for(int j = n_nodal_types.size() - 2; j > 0; --j){
      reps.push_back(std::accumulate(std::begin(n_nodal_types), std::begin(n_nodal_types) + j, 1, std::multiplies<int>()));
    }
    std::reverse(reps.begin(), reps.end());

    //component nodal types of causal type i
    std::vector<std::string> ct_i;
    ct_i.reserve(n_nodes);

    for(int j = 0; j < n_nodes; ++j){
      int pos = std::ceil((i+1)/double(reps[j]));
      pos = pos % n_nodal_types[j+1];

      if(pos == 0){
        pos = n_nodal_types[j+1];
      }

      std::vector<std::string> nt_j = nodal_types[j];
      ct_i.push_back(nt_j[pos - 1]);
    }
    //apply
  }
}
