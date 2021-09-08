/*
 * Functions for calculating genomic estimates of 
 * mean cross performance based on a GS model fitting 
 * additive and dominance (digenic) effects. Effects 
 * are modeled using AlphaSimR's coding scheme. 
 */

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

/*
 * Calculate cross means for different ploidy levels 
 * Input: 
 *  geno, genotypes (0,1,..,ploidy) 
 *        rows for indiviuals 
 *        columns for SNPs 
 *  a, additive effects for each SNP 
 *  d, dominance effects for each SNP 
 *  ploidy, can be 2, 4, or 6 
 * Output: 
 *  A matrix for all half-diallel crosses 
 *    Columns 1 and 2 are the parents of the cross 
 *    Column 3 is the expected mean of the cross 
 */
// [[Rcpp::export]]
arma::mat calcCrossMean(arma::mat& geno,
                        arma::vec& a,
                        arma::vec& d,
                        arma::uword ploidy){
  
  // Determine numbers and allocate output matrix
  arma::uword nInd = geno.n_rows;
  arma::uword nSnp = geno.n_cols;
  arma::uword nCrosses = nInd * (nInd-1) / 2;
  arma::mat output(nCrosses, 3, arma::fill::zeros);
  
  // Convert genotype data to a more convient type and layout
  arma::umat genoT = arma::conv_to<arma::umat>::from(geno.t());
  
  // Create additive and dominance genotype effects vectors
  // Uses AlphaSimR coding scheme
  arma::vec x = arma::regspace(0, ploidy); // Raw genotype
  x /= double(ploidy); // Proportional genotype
  arma::vec xa = 2.0*x - 1.0; // Additive effect
  arma::vec xd = -4.0*(x%x) + 4.0*x; // Digenic dominance effect
  
  // Calculate gamete probabilities
  // Probabilities determined by considering all possible combinations
  arma::mat gam;
  if(ploidy==2){
    gam = {
      {2, 0}, // 0 genotype (0, 1 gametes)
      {1, 1}, // 1
      {0, 2}  // 2
    };
    gam /= 2.0; // 2 choose 1
  }else if(ploidy==4){
    /*
     * A matrix of tetraploid gamete probabilities 
     * Assumes independent assortment of chromosomes
     * e.g. bivalent pairing only
     * The probabilities don't change too much with 
     * quadrivalents, so these seem to be reasonable
     * values even when the assumption of independent 
     * assortment of chromosomes is violated.
     */
    gam = {
      {6, 0, 0}, // 0 genotype (0, 1, 2 gametes)
      {3, 3, 0}, // 1 
      {1, 4, 1}, // 2 
      {0, 3, 3}, // 3 
      {0, 0, 6}  // 4 
    };
    gam /= 6.0; // 4 choose 2
  }else if(ploidy==6){
    gam = {
      {20,  0,  0,  0}, // 0 genotype (0, 1, 2, 3 gametes)
      {10, 10,  0,  0}, // 1 
      { 4, 12,  4,  0}, // 2 
      { 1,  9,  9,  1}, // 3 
      { 0,  4, 12,  4}, // 4 
      { 0,  0, 10, 10}, // 5 
      { 0,  0,  0, 20}  // 6 
    };
    gam /= 20.0; // 6 choose 3
  }else{
    Rcpp::stop("No gamete probabilities for this ploidy");
  }
  
  // Map parental genotypes to mean additive and dominance effects of progeny
  arma::mat mapA(ploidy+1, ploidy+1, arma::fill::zeros);
  arma::mat mapD(ploidy+1, ploidy+1, arma::fill::zeros);
  
  // Loop over parental genotypes to
  // compute means for each pair
  for(arma::uword i=0; i<=ploidy; i++){
    for(arma::uword j=0; j<=ploidy; j++){
      
      // Determine frequency of gamete pairs
      arma::mat F = gam.row(i).t() * gam.row(j);
      
      // Loop over gamete pairs
      for(arma::uword k=0; k<=(ploidy/2); k++){
        for(arma::uword l=0; l<=(ploidy/2); l++){
          
          // Multiply frequency by effect
          mapA(i,j) += F(k,l) * xa(k+l);
          mapD(i,j) += F(k,l) * xd(k+l);
        }
      }
      
    }
  }
  
  // Calculate cross means for partial diallel
  arma::uword k=0; // Cross identifier
  for(arma::uword i=0; i<(nInd-1); i++){
    for(arma::uword j=i+1; j<nInd; j++){
      
      // Record parents
      output(k,0) = i;
      output(k,1) = j;
      
      // Sum effects for all loci
      for(arma::uword m=0; m<nSnp; m++){
        output(k,2) +=
          mapA(genoT(m,i), genoT(m,j)) * a(m) +
          mapD(genoT(m,i), genoT(m,j)) * d(m);
      }
      
      k++;
    }
  }
  
  // C++ to R indexing
  output.col(0) += 1.0; 
  output.col(1) += 1.0; 
  
  return output;
}

