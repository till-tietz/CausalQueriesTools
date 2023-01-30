#include <RcppArmadillo.h>
#include <numeric>
#include <math.h>
#include <unistd.h>
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
void query_to_ct_c(std::vector<std::string> nodes,
                   std::vector<std::string> pars, //paste0(model$parameters_df$node, model$parameters_df$nodal_type) >> ensures correct handling under confounding
                   List nodal_types,
                   int n_causal_types,
                   int n_cores){

  int n_nodes = nodes.size();

  //data structures for component nodal types of a causal type -----------------
  std::vector<int> n_nodal_types(n_nodes);

  for(int i = 0; i < n_nodes; ++i) {
    std::vector<std::string> nodal_types_i = nodal_types[i];
    n_nodal_types[i] = nodal_types_i.size();
  }

  n_nodal_types.insert(n_nodal_types.begin(),1);
  n_nodal_types.insert(n_nodal_types.end(),1);
  std::vector<int> reps;

  for(int i = n_nodal_types.size() - 2; i > 0; --i) {
    reps.push_back(std::accumulate(std::begin(n_nodal_types), std::begin(n_nodal_types) + i,
                                   1, std::multiplies<int>()));
  }
  std::reverse(reps.begin(), reps.end());

  //loop over causal types in parallel -----------------------------------------
  #if defined(_OPENMP)
    #pragma omp parallel num_threads(n_cores)
    #pragma omp for
  #endif

  for(int i = 0; i < n_causal_types; i++){

    //generate causal type components ------------------------------------------
    //component nodal types of causal type i
    std::vector<std::string> ct_i;
    ct_i.reserve(n_nodes);
    //component nodal types with node name (to generate P matrix)
    std::vector<std::string> ct_i_named;
    ct_i_named.reserve(n_nodes);

    for(int j = 0; j < n_nodes; ++j){
      int pos = std::ceil((i+1)/double(reps[j]));
      pos = pos % n_nodal_types[j+1];

      if(pos == 0){
        pos = n_nodal_types[j+1];
      }

      std::vector<std::string> nt_j = nodal_types[j];
      ct_i.push_back(nt_j[pos - 1]);
      ct_i_named.push_back(nodes[j] + nt_j[pos - 1]);
    }

    //get P for causal type (map from parameters into causal types) ------------
    std::vector<int> P(pars.size());
    for(int j = 0; j < pars.size(); ++j) {
      auto pos = std::find(ct_i_named.begin(), ct_i_named.end(), pars[i]);
      if(pos != ct_i_named.end()) {
        P[i] = 1;
      }
    }

  }
}
