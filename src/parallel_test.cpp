#include <unistd.h>
#include <Rcpp.h>

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
int add_par(int n, int ncores) {
  int res = 0;

  #if defined(_OPENMP)
    #pragma omp parallel num_threads(ncores)
    #pragma omp for
  #endif

  for(int i = 0; i < n; i++) {
    sleep(1);
    res += 1;
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


