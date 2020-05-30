#include <Rcpp.h>
#include <random>
#include <omp.h>
#include "ABCoptim_types.h"

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

/**
 * This temp structure is used to save memory allocation
 */

#define UNIF_PARALLEL(a) runif_funs[(a)]->operator()((*rand_engines[(a)]))

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
  
  /**
   * OpenMP variables
   */
  unsigned int ncores;
  unsigned int seed;
  std::vector< int > pllblock;
  std::vector< std::mt19937_64 * > rand_engines;
  std::vector< std::uniform_real_distribution< double > * > runif_funs;
  
  // Setting up the RNG
  // - The first line creates an engine that uses the 64-bit Mersenne Twister by
  //   Matsumoto and Nishimura, 2000. One seed per core.
  // - The second line creates a function based on the real uniform between -1
  //   and 1. This receives as argument the engine
  // std::mt19937_64 engine((core_num + seed)*10);
  // std::uniform_real_distribution<double> my_runif(0.0, 1.0);

  /**Aux variables
   * 
   */
  unsigned int counts = 0u;
  double GlobalMinOld;
  std::vector< unsigned int > param2change;
  unsigned int neighbour = 0u;
  std::vector< double > ObjValSol;
  double FitnessSol = 0.0;
  // unsigned int t = 0;
  // unsigned int i = 0;
  std::vector< Input_type > solution;
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
    unsigned int criter_,
    unsigned int ncores_,
    unsigned int seed_
    );
  ~BeeHive() {
    for (unsigned int i = 0u; i < ncores; ++i) {
      delete rand_engines[i];
      delete runif_funs[i];
    }
  };
  
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
  unsigned int criter_,
  unsigned int ncores_,
  unsigned int seed_
) : fn(fn_), lb(lb_), ub(ub_), FoodNumber(FoodNumber_),
  Npar(GlobalParams_.size()), limit(limit_),
  maxCycle(maxCycle_), criter(criter_),
  Foods(FoodNumber_), f(FoodNumber_), prob(FoodNumber_), trial(FoodNumber_),
  fitness(FoodNumber_), GlobalParams(GlobalParams_), solution(FoodNumber_, GlobalParams_) {
  
  for (unsigned int i = 0u; i < FoodNumber; ++i)
    Foods[i] = clone(GlobalParams_);
  
  // Making room
  ans.reserve(maxCycle);
  
  /**
   * Preparing parallel computing
   */
  ncores = ncores_;
  seed   = seed_;
  rand_engines.reserve(ncores);
  runif_funs.reserve(ncores);
  ObjValSol.resize(FoodNumber, 0.0);
  param2change.resize(FoodNumber, 0u);
  
  pllblock.resize(ncores + 1, floor(FoodNumber / ncores));
  
  // Figuring out the size of the blocks
  pllblock[0u] = 0;
  pllblock[1] += (FoodNumber - std::accumulate(pllblock.begin(), pllblock.end(), 0));
  
  for (unsigned int i = 0u; i < ncores; ++i) {
    rand_engines.push_back(new std::mt19937_64(i + seed));
    runif_funs.push_back(new std::uniform_real_distribution<double>(0.0, 1.0));
    // std::cout << "Engine number " << i << " done." << std::endl;
    pllblock[i + 1] += pllblock[i];
    // printf("Block %i iterates in %i->%i\n", i, pllblock[i], pllblock[i+1u]);
  }
  
  omp_set_num_threads(ncores);
  
  // Initializing
  GlobalMin = as<double>(fn(GlobalParams));
  
  // Should be distributed equally
  for (unsigned int i = 0u; i < FoodNumber; ++i) {
    
    for (unsigned int j = 0u; j < Npar; ++j) {
      Foods[i][j] = lb[j] + (ub[j] - lb[j])/((double) FoodNumber - 1.0) * ((double) i);
      // Rprintf("Iteration Foods(%i, %j)=%.4f\n", i, j,Foods[i][j]);
    }
    
    // Checking if it is defined
    f[i] = as<double>(fn(Foods[i]));
    if (std::isnan(f[i]))
      throw std::logic_error("Undefined value for -fn-. Check the function's domain.");
    
    fitness[i] = CalculateFitness(f[i]);
    trial[i]   = 0;
  } 
  
  // for (auto i = Foods.begin(); i!=Foods.end();++i)
  //   print(*i);
  // stop("ASda");
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

//   if (Npar > 10u) {
//     
//     int core_num = -1;
// #pragma omp parallel for shared(Foods, ub, lb, index) default(none) \
//     firstprivate(core_num, Npar)
//     for (unsigned int j = 0u; j < Npar; ++j) {
//       if (core_num < 0)
//         core_num = omp_get_thread_num();
//       
//       Foods[index][j] = UNIF_PARALLEL(core_num) * (ub[j] - lb[j]) + lb[j];
//     }
//     
//   } else {
    
    for (unsigned int j = 0u; j < Npar; ++j)
      Foods[index][j] = UNIF_PARALLEL(0u) * (ub[j] - lb[j]) + lb[j];
    
  // }
  
  
  f[index]       = as<double>(fn(Foods[index]));
  fitness[index] = CalculateFitness(f[index]);
  trial[index]   = 0;
  
  return;
}

// Sends employed Bees
template <typename Input_type>
inline void BeeHive<Input_type>::SendEmployedBees() {
  
  // In principle, this function can be optimized
  int core_num = -1;
  double * solutionptr;
  double * Foodsptr;
  double * neighbourptr;
  
#pragma omp parallel private(neighbour, neighbourptr, FitnessSol, solutionptr, Foodsptr) \
  default(none) firstprivate(core_num, FoodNumber, Npar) \
  shared(Foods, lb, ub, ObjValSol, pllblock, solution, param2change)
  {
    
    core_num = omp_get_thread_num();

    // Generating the random proposals
    for (unsigned int i = pllblock[core_num]; i < pllblock[core_num + 1]; ++i) {
      
      solutionptr = (double *) &solution[i][0u];
      Foodsptr    = (double *) &Foods[i][0u];
  
      // Random parameter to change
      param2change[i] = (int) (UNIF_PARALLEL(core_num) * Npar);
      
      // Random neighbour to select
      neighbour = (int) (UNIF_PARALLEL(core_num) * (FoodNumber - 1u));;
      if (neighbour >= i)
        neighbour++;
      
      neighbourptr = (double *) &Foods[neighbour][0u];
      

      // Suggesting new solution
      for (unsigned int ii = 0; ii < Npar; ++ii)
        solutionptr[ii] = Foodsptr[ii];
      
      solutionptr[param2change[i]] = Foodsptr[param2change[i]] + 
        (
            Foodsptr[param2change[i]] -
              neighbourptr[param2change[i]]
        ) * (UNIF_PARALLEL(core_num) - .5) * 2.0;
      
      // Truncating
      if (solutionptr[param2change[i]] < lb[param2change[i]])
        solutionptr[param2change[i]] = lb[param2change[i]];
      
      if (solutionptr[param2change[i]] > ub[param2change[i]])
        solutionptr[param2change[i]] = ub[param2change[i]];
      
    }

#pragma omp barrier
    {
      if (core_num == 0u) {
        for (unsigned int i = 0u; i < FoodNumber; ++i) {
          ObjValSol[i]  = as<double>(fn(solution[i]));
          // printf("ObjValSol[%i]= %.4f, [", i, ObjValSol[i]);
          // for (int ii = 0; ii < Npar; ++ii)
          //   printf("%.4f, ", solution[i][ii]);
          // printf("\n");
        }
      }
    }
#pragma omp barrier

    for (unsigned int i = pllblock[core_num]; i < pllblock[core_num + 1]; ++i) {
      
      FitnessSol = CalculateFitness(ObjValSol[i]);
  
      if (FitnessSol > fitness[i]) {
        
        solutionptr = (double *) &solution[i][0u];
        Foodsptr    = (double *) &Foods[i][0u];
        
        Foodsptr[param2change[i]] = solutionptr[param2change[i]];
        // std::swap(Foods[i], solution[core_num]);
        fitness[i] = FitnessSol;
        f[i]       = ObjValSol[i];
        trial[i]   = 0u;
        
      } else 
        ++trial[i];
    }
  }
   
  return;
}

template <typename Input_type>
inline void BeeHive<Input_type>::CalculateProbabilities() {
  maxfit = fitness[0u];
  for (unsigned int i = 1u; i < FoodNumber; ++i) 
    if (fitness[i] > maxfit)
      maxfit = fitness[i];
   
#pragma omp parallel for simd
  for (unsigned int i = 0u; i < FoodNumber; ++i) 
    prob[i] = (0.9*((fitness[i])/(maxfit + 1e-20))) + 0.1;
  
  return;
  
}

template <typename Input_type>
inline void BeeHive<Input_type>::SendOnlookerBees() {
  
  unsigned int t = 0u, i = 0u;
  while (t < FoodNumber) {
    // Randomly choose a food source
    if (UNIF_PARALLEL(0u) < prob[i]) {
      t++;
      
      // Random parameter to change
      param2change[0u] = (unsigned int)(UNIF_PARALLEL(0u) * Npar);
      
      // Random neighbour to select
      neighbour = (unsigned int)(UNIF_PARALLEL(0u)*(FoodNumber - 1u));
      if (neighbour >= i)
        neighbour++;
      
      // Suggesting new solution
      solution[0u] = clone(Foods[i]);
      solution[0u][param2change[0u]] = Foods[i][param2change[0u]] + 
        (
            Foods[i][param2change[0u]] -
              Foods[neighbour][param2change[0u]]
        )*(UNIF_PARALLEL(0u) - 0.5) * 2.0;
      
      // Truncating
      if (solution[0u][param2change[0u]] < lb[param2change[0u]])
        solution[0u][param2change[0u]] = lb[param2change[0u]];
      
      if (solution[0u][param2change[0u]] > ub[param2change[0u]])
        solution[0u][param2change[0u]] = ub[param2change[0u]];
      
      // Comparing current solution with new one
      ObjValSol[0u]  = as<double>(fn(solution[0u]));
      FitnessSol = CalculateFitness(ObjValSol[0u]);
      
      if (FitnessSol > fitness[i]) {
        
        Foods[i][param2change[0u]] = solution[0u][param2change[0u]];
        fitness[i] = FitnessSol;
        f[i]       = ObjValSol[0u];
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

// [[Rcpp::export(name=".abc_cpp", rng = false)]]
List abc_cpp(
    NumericVector & par, 
    Function & fn,
    const std::vector< double > & lb,
    const std::vector< double > & ub,
    int FoodNumber   = 20, 
    int limit        = 100,
    int maxCycle     = 1000,
    int criter       = 50,
    int ncores       = 2,
    int seed         = 1 // , double tol=1e-10
) {
  
  // Initialize:
  BeeHive<NumericVector> bees(fn, par, lb, ub, FoodNumber, limit, maxCycle, criter, ncores, seed);

// #pragma omp parallel for shared(bees) firstprivate(ncores)
//   for (unsigned int i = 0u; i < ncores; ++i) {
//     double ans = bees.runif_funs[omp_get_thread_num()]->operator()((*bees.rand_engines[omp_get_thread_num()]));
//     std::cout << "Core number " << i << " generated " << ans << std::endl;
//   }
//   
//   return List::create();
  
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

ans <- .abc_cpp(c(1,1),fun, ub = c(5,5),lb = c(-5,-5), criter = 100, maxCycle = 100, ncores = 3)

stop()

library(microbenchmark)
library(ABCoptim)
microbenchmark(
  .abc_cpp(c(1,1),fun, ub = c(5,5),lb = c(-5,-5), criter = 100, maxCycle = 100, ncores = 3, FoodNumber = 100),
  abc_optim(c(1,1), fun, ub = 5, lb=-5, criter=100, maxCycle = 100, FoodNumber = 100),
  # abc_cpp(c(1,1), fun, ub = 5, lb=-5, criter=100, maxCycle = 100, FoodNumber = 20),
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