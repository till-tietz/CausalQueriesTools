#include <unistd.h>
#include <vector>
#include <omp.h>
#include <iostream>
#include <chrono>
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>

// [[Rcpp::export]]
int add_par(int n, int ncores) {
  int res = 0;
  #pragma omp parallel num_threads(ncores) reduction(+: res)
  {
    #pragma omp for
    for(int i = 0; i < n; i++) {
      sleep(1);
      res += 1;
    }
  }
  return res;
}


// [[Rcpp::export]]
int add_seq(int n){
  int res = 0;

  for(int i = 0; i < n; i++) {
    sleep(1);
    res += 1;
  }

  return res;
}



// [[Rcpp::export]]
std::vector<int> write_vec_par(std::vector<int> x, int y, int ncores) {
  std::vector<int> res(x.size());

  // Copy elements of x into res
  std::copy(x.begin(), x.end(), res.begin());
  Progress p(x.size(), true);

  #pragma omp parallel num_threads(ncores) default(none) shared(res,y,p)
  {
    #pragma omp for
    for(int i = 0; i < res.size(); i++) {
      if ( ! Progress::check_abort() ) {
        p.increment();
        res[i] = res[i] + y;
      }
    }
  }
  return res;
}

