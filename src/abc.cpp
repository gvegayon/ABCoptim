#include <Rcpp.h>
// #include <omp.h>
#include "ABCoptim_types.h"

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

/***
 * This temp structure is used to save memory allocation
 */


template <typename Input_type>
class BeeHive {
public:
  
  /**
   * Main arguments
   */
  Function & fn;
  const std::vector< double > & lb;
  const std::vector< double > & ub;
  
  /**Optim parameters
   * 
   */
  const unsigned int FoodNumber;
  const unsigned int Npar;
  const unsigned int limit;
  const unsigned int maxCycle;
  const unsigned int criter;
  
  /**Food sources
   * 
   */
  std::vector< Input_type > Foods;
  std::vector< double > f;
  std::vector< double > prob;
  std::vector< unsigned int > trial;
  std::vector< double > fitness;
  
  /**@brief Current optimum parameters
   * 
   */
  Input_type GlobalParams;
  double GlobalMin;
  unsigned int unchanged;
  std::vector< Input_type > ans;

  /**Aux variables
   * 
   */
  unsigned int counts = 0u;
  double GlobalMinOld;
  unsigned int param2change = 0u;
  unsigned int neighbour = 0u;
  double ObjValSol = 0.0;
  double FitnessSol = 0.0;
  // unsigned int t = 0;
  // unsigned int i = 0;
  Input_type solution;
  unsigned int maxtrialindex = 0u;
  double maxfit;
  
public:
  
  /**
   * Constructor and destructor
   */
  BeeHive(
    Function & fn_,
    const Input_type & GlobalParams_,
    const std::vector< double > & lb_,
    const std::vector< double > & ub_,
    unsigned int FoodNumber_,
    unsigned int limit_,
    unsigned int maxCycle_,
    unsigned int criter_
    );
  ~BeeHive() {};
  
  /**
   * Upate functions
   */
  
  
  double CalculateFitness(const double & x);
  void MemorizeBetsSource();
  void init(const unsigned int & index);
  void SendEmployedBees();
  void CalculateProbabilities();
  void SendOnlookerBees();
  void SendScoutBees();
  void run();
  
};
 
template <typename Input_type>
inline BeeHive<Input_type>::BeeHive(
  Function & fn_,
  const Input_type & GlobalParams_,
  const std::vector< double > & lb_,
  const std::vector< double > & ub_,
  unsigned int FoodNumber_,
  unsigned int limit_,
  unsigned int maxCycle_,
  unsigned int criter_
) : fn(fn_), lb(lb_), ub(ub_), FoodNumber(FoodNumber_),
  Npar(GlobalParams_.size()), limit(limit_),
  maxCycle(maxCycle_), criter(criter_),
  Foods(FoodNumber_, GlobalParams_), f(FoodNumber_), prob(FoodNumber_), trial(FoodNumber_),
  fitness(FoodNumber_), GlobalParams(GlobalParams_) {
  
  // Making room
  ans.reserve(maxCycle);
  
  // Initializing
  GlobalMin = as<double>(fn(GlobalParams));
  
  // Should be distributed equally
  for (unsigned int i = 0u; i < FoodNumber; ++i) {
    // Rprintf("Iteration %i\n", i);
    
    for (unsigned int j = 0u; j < Npar; ++j)
      Foods[i][j] = lb[j] + (ub[j] - lb[j])/(FoodNumber - 1.0)*i;
    
    // Checking if it is defined
    f[i] = as<double>(fn(Foods[i]));
    if (std::isnan(f[i]))
      throw std::logic_error("Undefined value for -fn-. Check the function's domain.");
    
    fitness[i] = CalculateFitness(f[i]);
    trial[i]   = 0;
  } 
  
  unchanged = 0;
  MemorizeBetsSource(); 
  ans.push_back(GlobalParams);
  unchanged = 0;
  
  return;
  
}

// Fitness computing
template <typename Input_type>
inline double BeeHive<Input_type>::CalculateFitness(const double & x) {
  return (x >= 0.0) ? 1.0/(x + 1.0) : 1.0 + std::abs(x);
}
 
// Sets the current best
template <typename Input_type>
inline void BeeHive<Input_type>::MemorizeBetsSource() {
  
  GlobalMinOld = GlobalMin;
  for (unsigned int i = 0u; i < FoodNumber; ++i)
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
template <typename Input_type>
inline void BeeHive<Input_type>::init(const unsigned int & index) {

  for (unsigned int j = 0u; j < Npar; ++j)
    Foods[index][j] = unif_rand() * (ub[j] - lb[j]) + lb[j];
  
  f[index]       = as<double>(fn(Foods[index]));
  fitness[index] = CalculateFitness(f[index]);
  trial[index]   = 0;
  
  return;
}

// Sends employed Bees
template <typename Input_type>
inline void BeeHive<Input_type>::SendEmployedBees() {
  
  // In principle, this function can be optimized
  for (unsigned int i = 0u; i < FoodNumber; ++i) {
    // Random parameter to change
    param2change = (int) (unif_rand() * Npar);
    
    // Random neighbour to select
    neighbour = (int) (unif_rand() * (FoodNumber - 1u));;
    if (neighbour >= i)
      neighbour++;
    
    // Suggesting new solution
    solution = clone(Foods[i]);
    solution[param2change] = Foods[i][param2change] + 
      (
          Foods[i][param2change] -
            Foods[neighbour][param2change]
    ) * (unif_rand() - .5) * 2;
    
    // Truncating
    if (solution[param2change] < lb[param2change])
      solution[param2change] = lb[param2change];
    
    if (solution[param2change] > ub[param2change])
      solution[param2change] = ub[param2change];
    
    // Comparing current solution with new one
    ObjValSol  = as<double>(fn(solution));
    FitnessSol = CalculateFitness(ObjValSol);
    if (FitnessSol > fitness[i]) {
      
      std::swap(Foods[i], solution);
      fitness[i] = FitnessSol;
      f[i]       = ObjValSol;
      trial[i]   = 0;
      
    } else 
      ++trial[i];
  } 
  
  return;
}

template <typename Input_type>
inline void BeeHive<Input_type>::CalculateProbabilities() {
  maxfit = fitness[0u];
  for (unsigned int i = 1u; i < FoodNumber; ++i) 
    if (fitness[i] > maxfit)
      maxfit = fitness[i];
   
// #pragma omp simd 
  for (unsigned int i = 0u; i < FoodNumber; ++i) 
    prob[i] = (0.9*((fitness[i])/(maxfit + 1e-20))) + 0.1;
  
  return;
  
}

template <typename Input_type>
inline void BeeHive<Input_type>::SendOnlookerBees() {
  
  unsigned int t = 0u, i = 0u;
  while (t < FoodNumber) {
    // Randomly choose a food source
    if (unif_rand() < prob[i]) {
      t++;
      
      // Random parameter to change
      param2change = (unsigned int)(unif_rand() * Npar);
      
      // Random neighbour to select
      neighbour = (unsigned int)(unif_rand()*(FoodNumber - 1u));
      if (neighbour >= i)
        neighbour++;
      
      // Suggesting new solution
      solution = clone(Foods[i]);
      solution[param2change] = Foods[i][param2change] + 
        (
            Foods[i][param2change] -
              Foods[neighbour][param2change]
        )*(unif_rand() - 0.5) * 2.0;
      
      // Truncating
      if (solution[param2change] < lb[param2change])
        solution[param2change] = lb[param2change];
      
      if (solution[param2change] > ub[param2change])
        solution[param2change] = ub[param2change];
      
      // Comparing current solution with new one
      ObjValSol  = as<double>(fn(solution));
      FitnessSol = CalculateFitness(ObjValSol);
      
      if (FitnessSol > fitness[i]) {
        
        std::swap(Foods[i], solution);
        fitness[i] = FitnessSol;
        f[i]       = ObjValSol;
        trial[i]   = 0u;
        
      } else 
        ++trial[i];
      
      
    } else { /* if */
      
      i++;
      if (i == FoodNumber)
        i=0u;
      
    }
  } /* while */
  
  return;
}

template <typename Input_type>
inline void BeeHive<Input_type>::SendScoutBees() {
  
  maxtrialindex = 0u;
  for(unsigned int i = 0u; i < FoodNumber; ++i) {
    if (trial[i] > trial[maxtrialindex])
      maxtrialindex = i;
  }
  
  // If it has reach the max, then init again
  if  (trial[maxtrialindex] >= limit)
    init(maxtrialindex);
  
  return;
}

template <typename Input_type>
inline void BeeHive<Input_type>::run() {

  counts = 0u;
  while (++counts < maxCycle) { 
    
    SendEmployedBees();
    CalculateProbabilities();
    SendOnlookerBees();
    MemorizeBetsSource();
    
    // Storing and breaking
    ans.push_back(GlobalParams);
    if (unchanged >= criter)
      break;
    
    SendScoutBees();
  }
  
  ans.shrink_to_fit();
  
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
  BeeHive<NumericVector> bees(fn, par, lb, ub, FoodNumber, limit, maxCycle, criter);

  // Optimize  
  bees.run();

  return List::create(
    _["Foods"]  = bees.Foods,
    _["f"]      = bees.f,
    _["fitness"]= bees.fitness,
    _["trial"]  = bees.trial,
    _["value"]  = bees.GlobalMin,
    _["par"]    = bees.GlobalParams,
    _["counts"] = bees.counts,
    _["hist"]   = bees.ans
  );
}


/***R
fun <- function(x) {
  -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2))
}

library(microbenchmark)
library(ABCoptim)
microbenchmark(
  # .abc_cpp(c(1,1),fun, ub = c(5,5),lb = c(-5,-5), criter = 100, maxCycle = 100)[-8],
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