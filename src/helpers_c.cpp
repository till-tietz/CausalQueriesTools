// [[Rcpp::depends(BH, bigmemory)]]
#include <RcppArmadillo.h>
#include <numeric>
#include <vector>
using namespace Rcpp;


//' Determines the number of times each nodes' vector of nodal types is repeated to generate the pattern of causal types.
//' Combining this information with the position of a causal type allows us to determine the nodal types used in the construction
//' of the given causal type.
//'
//' parameter nodal_types a List of nodal types
//' parameter nodes string vector of node names
//' return a map of node names and an integer of repetitions for each nodes' nodal type vector

std::map<std::string, int> get_causal_type_pattern(List &nodal_types,
                                                   std::vector<std::string> &nodes){

  //number of nodal types per node
  std::vector<int> n_nodal_types(nodal_types.size());

  for(int i = 0; i < nodal_types.size(); ++i){
    std::vector<std::string> nodal_types_i = nodal_types[i];
    n_nodal_types[i] = nodal_types_i.size();
  }

  n_nodal_types.insert(n_nodal_types.begin(),1);
  n_nodal_types.insert(n_nodal_types.end(),1);
  std::vector<int> reps;

  for(int i = n_nodal_types.size() - 2; i > 0; --i){
    reps.push_back(std::accumulate(std::begin(n_nodal_types), std::begin(n_nodal_types) + i, 1, std::multiplies<int>()));
  }
  std::reverse(reps.begin(), reps.end());


  std::map<std::string, int> map;
  std::transform(nodes.begin(), nodes.end(), reps.begin(), std::inserter(map, map.end()),
                 [](std::string const &s, int i) {
                   return std::make_pair(s, i);
                 });

  return map;
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




