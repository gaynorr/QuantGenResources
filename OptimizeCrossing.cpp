// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <omp.h>

// Find row given mapping index
// k = mapping index
// n = dimension of matrix
arma::uword mapRow(const arma::uword& k, 
                   const arma::uword& n){
  return n-2-static_cast<arma::uword>(sqrt(-8*double(k) + 4*double(n)*(double(n)-1)-7)/2-0.5);
}

// Find column given mapping index and row
// row = row from mapRow
// k = mapping index
// n = dimension of matrix
arma::uword mapCol(const arma::uword& row, 
                   const arma::uword& k, 
                   const arma::uword& n){
  return k+row+1 - n*(n-1)/2 + (n-row)*((n-row)-1)/2;
}

// Randomly samples integers without replacement
// nCross = number of integers to return
// nInd = number of integers to sample from
// Returns an integer vector of length nCross with values ranging from 0 to nInd-1
// Uses Jeffrey Scott Vitter's Method D
arma::uvec sampleInt(arma::uword n, arma::uword N){
  arma::uvec output;
  output.set_size(n);
  if(n == 0){
    return output;
  }
  double q, v, x, y1, y2;
  arma::uword threshold = 13*n;
  arma::uword S, limit, top, bottom;
  arma::vec u(1,arma::fill::randu);
  v = exp(log(u(0))/double(n));
  q = double(N-n+1);
  while((n>1) & (threshold<N)){
    while(true){
      while(true){
        x = double(N)*(1-v);
        S = floor(x);
        if(double(S)<q){
          break;
        }
        u.randu();
        v = exp(log(u(0))/double(n));
      }
      u.randu();
      y1 = exp(log(u(0)*double(N)/q)/double(n-1));
      v = y1*(1-x/double(N))*(q/(q-double(S)));
      if(v <= 1){
        break;
      }
      y2 = 1;
      top = N-1;
      if((n-1) > S){
        bottom = N-n;
        limit = N-S;
      }else{
        bottom = N-S-1;
        limit = N-n+1;
      }
      for(arma::uword i=N-1; i>=limit; --i)
        y2 *= double(top)/double(bottom);
      u.randu();
      if((double(N)/(double(N)-x)) >= (y1*exp(log(y2)/double(n-1)))){
        v = exp(log(u(0))/double(n-1));
        break;
      }
      v = exp(log(u(0))/double(n));
    }
    output(n-1) = S+1;
    N = N-S-1;
    --n;
    q = double(N-n+1);
    threshold -= 13;
  }
  if(n > 1){
    top = N-n;
    while(n >= 2){
      u.randu();
      S = 0;
      q = double(top)/double(N);
      while(q > u(0)){
        ++S;
        --top;
        --N;
        q = (q*double(top))/double(N);
      }
      output(n-1) = S+1;
      --N;
      --n;
    }
    u.randu();
    output(0) = floor(u(0)*N);
  }else{
    output(0) = floor(v*N);
  }
  return cumsum(output);
}

// Samples half diallele combinations without replacement
// If more than the total number of combinations is requested,
// only the number of combinations beyond the complete number is sampled
// nLevel = number of parents
// n = number of combinations to sample
arma::umat sampHalfDialComb(arma::uword nLevel, 
                            arma::uword n){
  arma::uword N = nLevel*(nLevel-1)/2;
  arma::uword fullComb = 0;
  // Determine number of complete combinations
  while(n>N){
    n -= N;
    ++fullComb;
  }
  arma::uvec samples = sampleInt(n,N);
  // Calculate selected combinations
  arma::umat output(2,n);
  for(arma::uword i=0; i<n; ++i){
    output(0,i) = mapRow(samples(i),nLevel);
    output(1,i) = mapCol(output(0,i),samples(i),nLevel);
  }
  // Add full combinations
  if(fullComb>0){
    arma::umat tmp(2,N*fullComb);
    arma::uword i;
    for(arma::uword j=0; j<(N*fullComb); ++j){
      i = j%N;
      tmp(0,j) = mapRow(i,nLevel);
      tmp(1,j) = mapCol(tmp(0,j),i,nLevel);
    }
    output = arma::join_rows(output,tmp);
  }
  return output;
}

// Finds parental contributions for a crossing plan
// crosses = the crossing plan
// nInd = number of potential parents
arma::vec calcContr(arma::uvec crosses, 
                    const arma::uword& nInd,
                    arma::uvec& rowPos, 
                    arma::uvec& colPos){
  double val = 1/(2*double(crosses.n_elem));
  arma::vec x(nInd, arma::fill::zeros);
  for(arma::uvec::iterator i=crosses.begin(); i!=crosses.end(); ++i){
    x(rowPos(*i)) += val;
    x(colPos(*i)) += val;
  }
  return x;
};

// Calculate the angle and length of a solution vector
// angle = angle of solution vector (output variable)
// length = length of solution vector (output variable)
// u = usefulness of solution
// sim = simularity of solution
// uMax = usefulness of maximum usefulness solution
// simMax = simularity of maximum usefulness solution
// uMin = usefulness of minimum simularity solution
// simMin = simularity of minimum simularity solution
void calcVec(double& angle, 
             double& length,
             double u, 
             double sim,
             const double& uMax, 
             const double& simMax,
             const double& uMin, 
             const double& simMin){
  u = (u-uMin)/(uMax-uMin);
  sim = (simMax-sim)/(simMax-simMin);
  length = sqrt(u*u+sim*sim);
  angle = acos(u/length);
  if(u<0){
    // Solution in wrong quadrant
    length = -length;
  }
}

// Performs mating between two crossing plan
// uses random selection of unique crosses
// a = a crossing plan
// b = a crossing plan (can have different length than a)
// returns crossing plan with same length as a
// Note: This implementation is slow and could be improved
arma::uvec mate(arma::uvec a, arma::uvec b){
  arma::uvec c = shuffle(unique(join_cols(a,b)));
  c.resize(a.n_elem);
  return c;
}

// Peforms mutation of a plan
// mutation is modeled as mating with a random crossing plan
// crosses = crossing plan to mutate (input/output variable)
// nMutate = length of random crossing plan
// potCross = number of potential crosses
void mutate(arma::uvec crosses, const arma::uword& nMutate, 
            const arma::uword& potCross){
  arma::uvec mutations = sampleInt(nMutate,potCross);
  crosses = mate(crosses,mutations);
}

// Quick and dirty implementation of a genetic algorithm (GA) 
// Uses usefulness criteria for crossing plan optimization
// Crosses are selected from a half diallele without selfs
// [[Rcpp::export]]
Rcpp::List selectCrosses(arma::uword nCross, //Number of crosses to make
                         double targetAngle, //Number of degrees to target from maximum usefulness (in radians)
                         arma::vec& u, //Vector of usefulness criterion for crosses
                         arma::mat& G, //Relationship matrix for individuals 
                         double probMut=0.01, //Mutation probability for progeny in GA
                         arma::uword nMutate=2, //Number of potential mutations in mutated progeny
                         arma::uword nSel=500, //Number of parents in GA
                         arma::uword nPop=10000, //Number of progeny in GA
                         arma::uword maxGen=1000, //Maximum number of generations
                         arma::uword maxRun=100, //Stopping criteria for maximum number of runs without change
                         double anglePenalty=0.5, //Penalty to vector length for off angle, higher value emphasises angle more
                         int nThreads=4){ //Number of threads for OpenMP
  omp_set_num_threads(nThreads); //Sets number of threads for OpenMP
  arma::umat crossPlan(2,nPop); //Crossing plan for GA generations
  arma::umat outCrossPlan(nCross,2); //Final crossing plan
  arma::uword nInd = G.n_cols; //Number of parents under consideration
  arma::uword potCross = nInd*(nInd-1)/2; //Number of potential crosses
  arma::umat Progeny(nCross,nPop), Parents(nCross,nSel); //Solutions
  arma::uvec Best(nCross); //Best solution
  arma::vec uProgeny(nPop), uParents(nSel); //Solution usefulness
  arma::vec simProgeny(nPop), simParents(nSel); //Solution simularity
  arma::vec angleProgeny(nPop), angleParents(nSel); //Solution angle
  arma::vec lenProgeny(nPop), lenParents(nSel); //Solution length
  arma::vec valProgeny(nPop), valParents(nSel); //Solution value
  arma::uvec rankProgeny(nPop); //Ranked progreny and select progeny
  double uBest, simBest, valBest, angleBest, lenBest; //Best solution usefulness, diversity, value, and angle
  double uMax, uMin, simMax, simMin; //Parameters for min simularity and max usefulness solutions
  arma::uword currentRun;
  arma::uvec rowPos(potCross), colPos(potCross); //Map cross number to row and column
  arma::uword tmpI=0;
  for(arma::uword i=0; i<(nInd-1); i++){
    for(arma::uword j=(i+1); j<nInd; j++){
      rowPos(tmpI) = i;
      colPos(tmpI) = j;
      ++tmpI;
    }
  }
  
  // Calculate maximum usefulness (doesn't require GA)
  arma::uvec uBestIndex = sort_index(u,"descend");
  uBestIndex.resize(nCross);
  arma::vec x = calcContr(uBestIndex.unsafe_col(0), nInd, 
                          rowPos, colPos);
  uMax = mean(u(uBestIndex));
  simMax = as_scalar(x.t()*G*x);
  
  if(targetAngle<1e-6){
    // No need to optimize crossing plan, just return best
    for(arma::uword i=0; i<nCross; ++i){
      outCrossPlan(i,0) = rowPos(uBestIndex(i));
      outCrossPlan(i,1) = colPos(uBestIndex(i));
    }
    return Rcpp::List::create(Rcpp::Named("crossPlan")=outCrossPlan+1, //C++ to R
                              Rcpp::Named("uMax")=uMax,
                              Rcpp::Named("simMax")=simMax);
  }
  
  // Calculate minimum simularity solution
  Rcpp::Rcout<<"Optimize for Simularity"<<std::endl<<std::endl;
  
  // Initialize progeny
#pragma omp parallel for schedule(static) 
  for(arma::uword i=0; i<nPop; i++){
    Progeny.col(i) = sampleInt(nCross, potCross);
    arma::vec x = calcContr(Progeny.unsafe_col(i), nInd, 
                            rowPos, colPos);
    simProgeny(i) = as_scalar(x.t()*G*x);
  }
  
  // Select parents and best solution
  rankProgeny = sort_index(simProgeny,"ascend");
  for(arma::uword i=0; i<nSel; i++){
    Parents.col(i) = Progeny.col(rankProgeny(i));
    simParents(i) = simProgeny(rankProgeny(i));
  }
  simBest = simParents(0);
  Best = Parents.col(0);
  
  // Run GA for simularity
  Rcpp::Rcout<<"Gen  Simularity"<<std::endl;
  currentRun = 0;
  for(arma::uword gen=0; gen<maxGen; gen++){
    // Mate parents
    crossPlan = sampHalfDialComb(nSel, nPop);
#pragma omp parallel for schedule(static)
    for(arma::uword i=0; i<nPop; i++){
      Progeny.col(i) = 
        mate(Parents.unsafe_col(crossPlan(0,i)),
             Parents.unsafe_col(crossPlan(1,i)));
      arma::vec p(1,arma::fill::randu);
      if(as_scalar(p)<probMut){
        mutate(Progeny.unsafe_col(i),nMutate,potCross);
      }
      arma::vec x = calcContr(Progeny.unsafe_col(i), nInd, 
                              rowPos, colPos);
      simProgeny(i) = as_scalar(x.t()*G*x);
    }
    // Select parents and best solution
    rankProgeny = sort_index(simProgeny,"ascend");
    for(arma::uword i=0; i<nSel; i++){
      Parents.col(i) = Progeny.col(rankProgeny(i));
      simParents(i) = simProgeny(rankProgeny(i));
    }
    if(simParents(0)<simBest){
      simBest = simParents(0);
      Best = Parents.col(0);
      currentRun = 0;
    }else{
      ++currentRun;
    }
    
    //Report status
    if(gen%10 == 0){
      Rcpp::Rcout<<gen<<"  "<<simBest<<std::endl;
    }
    
    //Terminate if solution isn't changing
    if(currentRun>=maxRun){
      break;
    }
  }
  simMin = simBest;
  uMin = mean(u(Best));
  
  // Perform optimization for best solution
  Rcpp::Rcout<<std::endl<<std::endl<<"Optimize for Crossing Plan"<<std::endl<<std::endl;
  
  // Initialize progeny
#pragma omp parallel for schedule(static) 
  for(arma::uword i=0; i<nPop; i++){
    Progeny.col(i) = sampleInt(nCross, potCross);
    arma::vec x = calcContr(Progeny.unsafe_col(i), nInd, 
                            rowPos, colPos);
    simProgeny(i) = as_scalar(x.t()*G*x);
    uProgeny(i) = mean(u(Progeny.col(i)));
    calcVec(angleProgeny(i), lenProgeny(i), uProgeny(i),
            simProgeny(i), uMax, simMax, uMin, simMin);
    valProgeny(i) = lenProgeny(i)-anglePenalty*fabs(angleProgeny(i)-targetAngle);
  }
  
  // Select parents and best solution
  rankProgeny = sort_index(valProgeny,"descend");
  for(arma::uword i=0; i<nSel; i++){
    Parents.col(i) = Progeny.col(rankProgeny(i));
    uParents(i) = uProgeny(rankProgeny(i));
    simParents(i) = simProgeny(rankProgeny(i));
    angleParents(i) = angleProgeny(rankProgeny(i));
    lenParents(i) = lenProgeny(rankProgeny(i));
    valParents(i) = valProgeny(rankProgeny(i));
  }
  valBest = valParents(0);
  angleBest = angleParents(0);
  lenBest = lenParents(0);
  uBest = uParents(0);
  simBest = simParents(0);
  Best = Parents.col(0);
  
  // Run GA
  Rcpp::Rcout<<"Gen  Usefulness  Simularity  Angle  Length  Value"<<std::endl;
  currentRun = 0;
  for(arma::uword gen=0; gen<maxGen; gen++){
    // Mate parents
    crossPlan = sampHalfDialComb(nSel, nPop);
#pragma omp parallel for schedule(static)
    for(arma::uword i=0; i<nPop; i++){
      Progeny.col(i) = 
        mate(Parents.unsafe_col(crossPlan(0,i)),
             Parents.unsafe_col(crossPlan(1,i)));
      arma::vec p(1,arma::fill::randu);
      if(as_scalar(p)<probMut){
        mutate(Progeny.unsafe_col(i),nMutate,potCross);
      }
      arma::vec x = calcContr(Progeny.unsafe_col(i), nInd, 
                              rowPos, colPos);
      simProgeny(i) = as_scalar(x.t()*G*x);
      uProgeny(i) = mean(u(Progeny.col(i)));
      calcVec(angleProgeny(i), lenProgeny(i), uProgeny(i),
              simProgeny(i), uMax, simMax, uMin, simMin);
      valProgeny(i) = lenProgeny(i)-anglePenalty*fabs(angleProgeny(i)-targetAngle);
    }
    
    // Select parents and best solution
    rankProgeny = sort_index(valProgeny,"descend");
    for(arma::uword i=0; i<nSel; i++){
      Parents.col(i) = Progeny.col(rankProgeny(i));
      uParents(i) = uProgeny(rankProgeny(i));
      simParents(i) = simProgeny(rankProgeny(i));
      angleParents(i) = angleProgeny(rankProgeny(i));
      lenParents(i) = lenProgeny(rankProgeny(i));
      valParents(i) = valProgeny(rankProgeny(i));
    }
    if(valParents(0)>valBest){
      valBest = valParents(0);
      angleBest = angleParents(0);
      lenBest = lenParents(0);
      uBest = uParents(0);
      simBest = simParents(0);
      Best = Parents.col(0);
      currentRun = 0;
    }else{
      ++currentRun;
    }
    
    //Report status
    if(gen%10 == 0){
      Rcpp::Rcout<<gen<<"  "<<uBest<<"  "<<simBest<<"  "<<angleBest<<"  "<<lenBest<<"  "<<valBest<<std::endl;
    }
    
    //Terminate if solution isn't changing
    if(currentRun>=maxRun){
      break;
    }
  }
  
  //Convert solution to an ordered crossing plan
  Best = sort(Best);
  for(arma::uword i=0; i<nCross; ++i){
    outCrossPlan(i,0) = rowPos(Best(i));
    outCrossPlan(i,1) = colPos(Best(i));
  }
  
  return Rcpp::List::create(Rcpp::Named("crossPlan")=outCrossPlan+1, //C++ to R
                            Rcpp::Named("uMax")=uMax,
                            Rcpp::Named("uMin")=uMin,
                            Rcpp::Named("simMax")=simMax,
                            Rcpp::Named("simMin")=simMin,
                            Rcpp::Named("uBest")=uBest,
                            Rcpp::Named("simBest")=simBest,
                            Rcpp::Named("angleBest")=angleBest,
                            Rcpp::Named("lenBest")=lenBest);
}
