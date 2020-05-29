#include <Rcpp.h>
#include <omp.h>
#include "ABCoptim_types.h"

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

/***
 * This temp structure is used to save memory allocation
 */


// Fitness computing
inline double CalculateFitness(const double & x) {
  return  (x >= 0.0) ? 1.0/(x + 1.0) : 1.0 + std::abs(x);
}

// Sets the current best
inline void MemorizeBetsSource(
    double & GlobalMin,
    NumericVector & GlobalParams,
    std::vector< double > & f, 
    std::vector< NumericVector > & Foods,
    int & unchanged
  ) {
  
  double GlobalMinOld = GlobalMin;
  for (unsigned int i = 0u; i < f.size(); ++i)
    if (f[i] < GlobalMin) {
      GlobalMin    = f[i];
      GlobalParams = clone(Foods[i]);
    }
  
  if (GlobalMin == GlobalMinOld)
    ++unchanged;
  else
    unchanged = 0;
  
  return;
}

// Initializes a food source
inline void init(
    int index,
    std::vector< double > & fitness,
    std::vector< double > & f,
    std::vector< int > & trial,
    Function & fun,
    std::vector< NumericVector > & Foods,
    const std::vector< double > & lb,
    const std::vector< double > & ub
  ) {

  for (unsigned int j = 0u; j < lb.size(); ++j)
    Foods[index][j] = unif_rand() * (ub[j] - lb[j]) + lb[j];
  
  f[index]       = as<double>(fun(Foods[index]));
  fitness[index] = CalculateFitness(f[index]);
  trial[index]   = 0;
  
  return;
}

// Sends employed Bees
inline void SendEmployedBees(
    double & GlobalMin,
    NumericVector & GlobalParams, 
    std::vector< double > & fitness,
    std::vector<double> & f,
    std::vector< int > & trial,
    Function & fun,
    std::vector< NumericVector > & Foods,
    const std::vector< double > & lb,
    const std::vector< double > & ub,
    TempParams & temp_pars
  ) {
  
  // NumericVector solution(Foods.ncol());
  for (unsigned int i = 0u; i < Foods.size(); ++i) {
    // Random parameter to change
    temp_pars.param2change = (int) (unif_rand() * temp_pars.solution.size());
    
    // Random neighbour to select
    temp_pars.neighbour = (int) (unif_rand() * (Foods.size() - 1u));;
    if (temp_pars.neighbour >= i)
      temp_pars.neighbour++;
    
    // Suggesting new solution
    temp_pars.solution = clone(Foods[i]);
    temp_pars.solution[temp_pars.param2change] = Foods[i][temp_pars.param2change] + 
      (
          Foods[i][temp_pars.param2change] -
            Foods[temp_pars.neighbour][temp_pars.param2change]
    ) * (unif_rand() - .5) * 2;
    
    // Truncating
    if (temp_pars.solution[temp_pars.param2change] < lb[temp_pars.param2change])
      temp_pars.solution[temp_pars.param2change] = lb[temp_pars.param2change];
    
    if (temp_pars.solution[temp_pars.param2change] > ub[temp_pars.param2change])
      temp_pars.solution[temp_pars.param2change] = ub[temp_pars.param2change];
    
    // Comparing current solution with new one
    temp_pars.ObjValSol  = as<double>(fun(temp_pars.solution));
    temp_pars.FitnessSol = CalculateFitness(temp_pars.ObjValSol);
    if (temp_pars.FitnessSol > fitness[i]) {
      
      std::swap(Foods[i], temp_pars.solution);
      fitness[i] = temp_pars.FitnessSol;
      f[i]       = temp_pars.ObjValSol;
      trial[i]   = 0;
      
    } else 
      ++trial[i];
  } 
  
  return;
}

inline void CalculateProbabilities(
    std::vector< NumericVector > & Foods,
    std::vector< double > & fitness,
    std::vector< double > & prob
  ) {
  double maxfit = fitness[0u];
  for (unsigned int i = 1u; i < Foods.size(); ++i) 
    if (fitness[i] > maxfit)
      maxfit = fitness[i];
  
#pragma omp simd 
  for (unsigned int i = 0u; i < Foods.size(); ++i) 
    prob[i] = (0.9*((fitness[i])/(maxfit + 1e-20))) + 0.1;
  
  return;
  
}


inline void SendOnlookerBees(
    double & GlobalMin,
    NumericVector & GlobalParams, 
    std::vector< double > & fitness,
    std::vector< double > & f,
    std::vector< int > & trial,
    std::vector< double > & prob,
    Function & fun,
    std::vector< NumericVector > & Foods,
    const std::vector< double > & lb,
    const std::vector< double > & ub,
    TempParams & temp_pars
  ) {
  
  unsigned int t = 0u, i = 0u;
  while (t < Foods.size()) {
    // Randomly choose a food source
    if (unif_rand() < prob[i]) {
      t++;
      
      // Random parameter to change
      temp_pars.param2change = (unsigned int)(unif_rand()*temp_pars.solution.size());
      
      // Random neighbour to select
      temp_pars.neighbour = (unsigned int)(unif_rand()*(Foods.size() - 1u));
      if (temp_pars.neighbour >= i)
        temp_pars.neighbour++;
      
      // Suggesting new solution
      temp_pars.solution = clone(Foods[i]);
      temp_pars.solution[temp_pars.param2change] = Foods[i][temp_pars.param2change] + 
        (
            Foods[i][temp_pars.param2change] -
              Foods[temp_pars.neighbour][temp_pars.param2change]
        )*(unif_rand() - 0.5) * 2.0;
      
      // Truncating
      if (temp_pars.solution[temp_pars.param2change] < lb[temp_pars.param2change])
        temp_pars.solution[temp_pars.param2change] = lb[temp_pars.param2change];
      
      if (temp_pars.solution[temp_pars.param2change] > ub[temp_pars.param2change])
        temp_pars.solution[temp_pars.param2change] = ub[temp_pars.param2change];
      
      // Comparing current solution with new one
      temp_pars.ObjValSol  = as<double>(fun(temp_pars.solution));
      temp_pars.FitnessSol = CalculateFitness(temp_pars.ObjValSol);
      
      if (temp_pars.FitnessSol > fitness[i]) {
        
        std::swap(Foods[i], temp_pars.solution);
        fitness[i] = temp_pars.FitnessSol;
        f[i]       = temp_pars.ObjValSol;
        trial[i]   = 0u;
        
      } else 
        ++trial[i];
      
      
    } else { /* if */
      
      i++;
      if (i == Foods.size())
        i=0u;
      
    }
  } /* while */
  
  return;
}

inline void SendScoutBees(
    std::vector< double > & fitness,
    std::vector< double > & f,
    std::vector< int > & trial,
    std::vector< double > & prob,
    Function & fun,
    std::vector< NumericVector > & Foods,
    const std::vector< double > & lb,
    const std::vector< double > & ub,
    int limit
) {
  
  unsigned int maxtrialindex = 0u;
  for(unsigned int i = 0u; i < Foods.size(); ++i) {
    if (trial[i] > trial[maxtrialindex])
      maxtrialindex = i;
  }
  
  // If it has reach the max, then init again
  if  (trial[maxtrialindex] >= limit)
    init(maxtrialindex, fitness, f, trial, fun, Foods, lb, ub);
  
  return;
}

// [[Rcpp::export(name=".abc_cpp")]]
List abc_cpp(
    NumericVector & par, 
    Function & fn,
    const std::vector< double > & lb,
    const std::vector< double > & ub,
    int FoodNumber   = 20, 
    int limit        = 100,
    int maxCycle     = 1000,
    int criter       = 50 // , double tol=1e-10
) {
  
  // Initialize:
  // NumericMatrix Foods(FoodNumber, par.size());
  std::vector< NumericVector > Foods(FoodNumber);
  NumericVector empty(par.size());
  for (auto iter = Foods.begin(); iter != Foods.end(); ++iter)
    *iter = clone(empty);
  
  // NumericVector Foods_row(par.size());
  
  std::vector< double > f(FoodNumber), fitness(FoodNumber), prob(FoodNumber);
  std::vector< int > trial(FoodNumber);
  
  NumericMatrix ans(maxCycle, par.size());

  // Initializing
  NumericVector GlobalParams = clone(par);
  double GlobalMin = as<double>(fn(GlobalParams));

  // Should be distributed equally
  for (int i = 0; i < FoodNumber; ++i) {
    // Rprintf("Iteration %i\n", i);
    
    for (int j = 0; j < par.size(); ++j)
      Foods[i][j] = lb[j] + (ub[j] - lb[j])/(FoodNumber - 1.0)*i;
    
    // Checking if it is defined
    f[i] = as<double>(fn(Foods[i]));
    if (NumericVector::is_na(f[i]))
      stop("Undefined value for -fn-. Check the function's domain.");
    
    fitness[i] = CalculateFitness(f[i]);
    trial[i]   = 0;
  } 
  
  int unchanged = 0;
  MemorizeBetsSource(GlobalMin, GlobalParams, f, Foods, unchanged); 
  ans(0,_) = clone(GlobalParams);
  unchanged = 0;

  // double OldGlobalMin = GlobalMin;
  int i=0;
  TempParams temp_pars(par.size());
  // NumericVector solution(par.size());
  while (++i<maxCycle) { 
    
    SendEmployedBees(GlobalMin, GlobalParams, fitness, f, trial, fn, Foods, lb, ub, temp_pars);
    // Rprintf("Ok employee\n");
    CalculateProbabilities(Foods, fitness, prob);
    // Rprintf("Ok prob\n");
    SendOnlookerBees(GlobalMin, GlobalParams, fitness, f, trial, prob, fn, Foods, lb, ub, temp_pars);
    // Rprintf("Ok onlooker\n");
    MemorizeBetsSource(GlobalMin, GlobalParams, f, Foods, unchanged);
    // Rprintf("Ok memorize\n"); 
    
    // Storing and breaking
    ans(i,_) = clone(GlobalParams);
    if (unchanged >= criter)
      break;
    
    // If not, then send them out again
    SendScoutBees(fitness, f, trial, prob, fn, Foods, lb, ub, limit);
  }
  
  return List::create(
    _["Foods"]  = Foods,
    _["f"]      = f,
    _["fitness"]= fitness,
    _["trial"]  = trial,
    _["value"]  = GlobalMin,
    _["par"]    = GlobalParams,
    _["counts"] = i,
    _["hist"]   = ans(Range(0,i),_)
  );
}


/***R
fun <- function(x) {
  -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2))
}

library(microbenchmark)
library(ABCoptim)
microbenchmark(
  .abc_cpp(c(1,1),fun, ub = c(5,5),lb = c(-5,-5), criter = 100, maxCycle = 100)[-8],
  abc_optim(c(1,1), fun, ub = 5, lb=-5, criter=100, maxCycle = 100, FoodNumber = 20),
  abc_cpp(c(1,1), fun, ub = 5, lb=-5, criter=100, maxCycle = 100, FoodNumber = 20),
  times=50,
  unit="relative"
)

# ans <- abc_cpp(c(1,1),fun, upper = 50,lower = -50, criter = 50, MNC = 500, SN = 20)
#
# fs <- function(x) sum(x^2)
# ans <- abc_cpp(rep(500,5), fs, lower = -10, upper=10, criter=200)
# tail(ans$ans)
#
# plot(ans$ans,type = "l")

*/