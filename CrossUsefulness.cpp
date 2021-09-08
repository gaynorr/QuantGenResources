// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <omp.h>

arma::uword mapRow(const arma::uword& k, const arma::uword& n){
  return n-2-static_cast<arma::uword>(sqrt(-8*double(k) + 4*double(n)*(double(n)-1)-7)/2-0.5);
}
arma::uword mapCol(const arma::uword& row, const arma::uword& k, const arma::uword& n){
  return k+row+1 - n*(n-1)/2 + (n-row)*((n-row)-1)/2;
}

// Calculates the usefulness of cross using simulated progeny
// [[Rcpp::export]]
arma::vec calcUsefulness(const arma::imat& parents, //Haplotypes (0-1)
                         const arma::imat& progeny, //Haplotypes (0-3)
                         const arma::vec& a,
                         int nBest=50,
                         int nThreads=4){
  omp_set_num_threads(nThreads);
  arma::uword nParents = parents.n_rows/2;
  arma::uword nProgeny = progeny.n_rows/2;
  arma::uword nLoci = parents.n_cols;
  arma::vec usefulness(nParents*(nParents-1)/2);
#pragma omp parallel for schedule(static) 
  for(arma::uword n=0; n<usefulness.n_elem; ++n){
    arma::uword i = mapRow(n,nParents);
    arma::uword j = mapCol(i,n,nParents);
    arma::vec blup(nProgeny, arma::fill::zeros);
    for(arma::uword k=0; k<nLoci; ++k){
      arma::vec x(4);
      x(0) = a(k)*parents(2*i,k);
      x(1) = a(k)*parents(2*i+1,k);
      x(2) = a(k)*parents(2*j,k);
      x(3) = a(k)*parents(2*j+1,k);
      for(arma::uword l=0; l<(2*nProgeny); ++l){
        blup(l/2) += x(progeny(l,k));
      }
    }
    blup = sort(blup,"descend");
    usefulness(n) = mean(blup(arma::span(0,nBest-1)));
  }
  return usefulness;
}
