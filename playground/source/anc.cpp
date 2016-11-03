#include <Rcpp.h>


using namespace Rcpp;

// Fitness computing
double CalculateFitness(double x) {
  return  (x>=0) ? 1/(x+1) : 1 + fabs(x);
}

// Sets the current best
void MemorizeBetsSource(
    double & GlobalMin, NumericVector & GlobalParams,
    NumericVector & f, NumericMatrix & Foods, int & unchanged) {
  
  int change=0;
  for (int i=0;i<f.size(); i++)
    if (f[i] < GlobalMin) {
      GlobalMin    = f[i];
      GlobalParams = Foods(i,_);
      ++change;
    }
  
  if (change>0) ++unchanged;
  
  return;
}

// Initializes a food source
void init(int index, NumericVector & fitness, NumericVector & f,
          IntegerVector & trial,
          Function & fun,
          NumericMatrix & Foods, double lb, double ub) {

  Foods(index,_) = runif(Foods.ncol(), lb, ub);
  
  f[index]       = as<double>(fun(Foods(index,_)));
  fitness[index] = CalculateFitness(f[index]);
  trial[index]   = 0;
  
  return;
}

// Sends employed Bees
void SendEmployedBees(
    double & GlobalMin, NumericVector & GlobalParams, 
    NumericVector & fitness, NumericVector & f,
    IntegerVector & trial,
    Function & fun,
    NumericMatrix & Foods, double lb, double ub) {
  
  int param2change, neighbour;
  double ObjValSol, FitnessSol;
  NumericVector solution(Foods.ncol());
  for (int i=0;i<Foods.nrow();i++) {
    // Random parameter to change
    param2change = (int)(unif_rand()*Foods.ncol());
    
    // Random neighbour to select
    neighbour    = i;
    while (neighbour == i) 
      neighbour = (int)(unif_rand()*Foods.nrow());
    
    // Suggesting new solution
    solution = Foods(i,_);
    solution[param2change] = Foods(i,param2change) + 
      (Foods(i, param2change) - Foods(neighbour, param2change))*(unif_rand()-.5)*2;
    
    // Truncating
    if (solution[param2change] < lb) solution[param2change] = lb;
    if (solution[param2change] > ub) solution[param2change] = ub;
    
    // Comparing current solution with new one
    ObjValSol  = as<double>(fun(solution));
    FitnessSol = CalculateFitness(ObjValSol);
    if (FitnessSol > fitness[i]) {
      Foods(i,_) = solution;
      fitness[i] = FitnessSol;
      f[i]       = ObjValSol;
      trial[i]   = 0;
    } else {
      trial[i]+=1;
    }
  } 
  
  return;
}

void CalculateProbabilities(NumericMatrix & Foods, NumericVector & fitness,
                            NumericVector & prob) {
  double maxfit = fitness[0];
  for (int i=0;i<Foods.nrow();i++) 
    if (fitness[i] > maxfit) maxfit = fitness[i];
    
  for (int i=0;i<Foods.nrow();i++) 
    prob[i] = (0.9*(fitness[i]/maxfit)) + 0.1;
  
  return;
  
}

void SendOnlookerBees(
    double & GlobalMin, NumericVector & GlobalParams, 
    NumericVector & fitness, NumericVector & f,
    IntegerVector & trial, NumericVector & prob,
    Function & fun,
    NumericMatrix & Foods, double lb, double ub) {
  
  int param2change, neighbour;
  double ObjValSol, FitnessSol;
  NumericVector solution(Foods.ncol());
  
  int t = 0, i=0;
  double r;
  while (t < Foods.nrow()) {
    // Randomly choose a food source
    // r = *Foods.nrow();
    if (unif_rand() < prob[i]) {
      t++;
      
      // Random parameter to change
      param2change = (int)(unif_rand()*Foods.ncol());
      
      // Random neighbour to select
      neighbour = i;
      while (neighbour == i) 
        neighbour = (int)(unif_rand()*Foods.nrow());
      
      // Suggesting new solution
      solution = Foods(i,_);
      solution[param2change] = Foods(i,param2change) + 
        (Foods(i, param2change) - Foods(neighbour, param2change))*(unif_rand()-.5)*2;
      
      // Truncating
      if (solution[param2change] < lb) solution[param2change] = lb;
      if (solution[param2change] > ub) solution[param2change] = ub;
      
      // Comparing current solution with new one
      ObjValSol  = as<double>(fun(solution));
      FitnessSol = CalculateFitness(ObjValSol);
      if (FitnessSol > fitness[i]) {
        Foods(i,_) = solution;
        fitness[i] = FitnessSol;
        f[i]       = ObjValSol;
        trial[i]   = 0;
      } else {
        trial[i]+=1;
      }
      
    } else { /* if */
      i++;
      if (i==Foods.nrow()) i=0;
    }
  } /* while */
  
  return;
}

void SendScoutBees(
    NumericVector & fitness, NumericVector & f,
    IntegerVector & trial, NumericVector & prob,
    Function & fun,
    NumericMatrix & Foods, double lb, double ub, int limit) {
  
  int maxtrialindex = 0;
  for(int i=1;i<Foods.nrow();i++) {
    if (trial[i] > trial[maxtrialindex])
      maxtrialindex = i;
  }
  
  // If it has reach the max, then init again
  if  (trial[maxtrialindex]>=limit)
    init(maxtrialindex, fitness, f, trial, fun, Foods, lb, ub);
  
  return;
}

// [[Rcpp::export]]
List abc_cpp(
    NumericVector & par, Function & fun,
    double upper= 1e15, double lower=-1e15,
    int SN = 20, int limit=100, int MNC=1000, int criter=50 // , double tol=1e-10
) {
  
  // Initialize:
  double ub = upper, lb=lower;

  NumericMatrix Foods(SN, par.size());        // Food sources
  NumericVector f(SN), fitness(SN), prob(SN); // Funval, Fitness Solutions amd Probs
  IntegerVector trial(SN);                    // Trial numbers
  
  NumericMatrix ans(MNC, par.size());

  // Initializing
  NumericVector GlobalParams = clone(par);
  double GlobalMin = as<double>(fun(GlobalParams));
  double num_eval = 0.0;
  for (int i=0;i<SN;i++) 
    init(i, fitness, f, trial, fun, Foods, lb, ub);

  int unchanged = 0;
  MemorizeBetsSource(GlobalMin, GlobalParams, f, Foods, unchanged);
  unchanged = 0;

  // Rprintf("%f",GlobalMin);
    
  // double OldGlobalMin = GlobalMin;
  int i=-1;
  while (++i<MNC) {
    
    SendEmployedBees(GlobalMin, GlobalParams, fitness, f, trial, fun, Foods, lb, ub);
    CalculateProbabilities(Foods, fitness, prob);
    SendOnlookerBees(GlobalMin, GlobalParams, fitness, f, trial, prob, fun, Foods, lb, ub);
    MemorizeBetsSource(GlobalMin, GlobalParams, f, Foods, unchanged);
    SendScoutBees(fitness, f, trial, prob, fun, Foods, lb, ub, limit);
    
    // Stop criteria
    ans(i,_) = GlobalParams;
    
    
    if (unchanged>=criter) break;
  }
  
  return List::create(
    _["Foods"]     = Foods,
    _["f"]         = f,
    _["fitness"]   = fitness,
    _["trial"]     = trial,
    _["GlobalMin"] = GlobalMin,
    _["GlobalParams"] = GlobalParams,
    _["counts"]    = i,
    _["ans"]       = ans
  );
}


/***R
fun <- function(x) {
  -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2))
}

library(microbenchmark)
library(ABCoptim)
microbenchmark(
  abc_cpp(c(1,1),fun, upper = 5,lower = -5, criter = 100, MNC = 100, SN = 20)[-8],
  abc_optim(c(1,1), fun, ub = 5, lb=-5, criter=100, maxCycle = 100, FoodNumber = 20), times=100,
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