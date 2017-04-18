#include <Rcpp.h>

using namespace Rcpp;

// Fitness computing
double CalculateFitness(double x) {
  return  (x>=0) ? 1/(x+1) : 1 + fabs(x);
}

// Sets the current best
void MemorizeBetsSource(
    double & GlobalMin,
    NumericVector & GlobalParams,
    NumericVector & f, 
    NumericMatrix & Foods,
    int & unchanged
  ) {
  
  double GlobalMinOld=GlobalMin;
  for (int i=0;i<f.size(); i++)
    if (f[i] < GlobalMin) {
      GlobalMin    = f[i];
      GlobalParams = Foods(i,_);
    }
  
  if (GlobalMin == GlobalMinOld) ++unchanged;
  else unchanged = 0;
  
  return;
}

// Initializes a food source
void init(
    int index,
    NumericVector & fitness,
    NumericVector & f,
    IntegerVector & trial,
    Function & fun,
    NumericMatrix & Foods,
    const NumericVector & lb,
    const NumericVector & ub
  ) {

  for (int j=0; j<lb.size(); j++)
    Foods.at(index,j) = runif(1, lb.at(j), ub.at(j))[0];
  
  f[index]       = as<double>(fun(Foods(index,_)));
  fitness[index] = CalculateFitness(f[index]);
  trial[index]   = 0;
  
  return;
}

// Sends employed Bees
void SendEmployedBees(
    double & GlobalMin,
    NumericVector & GlobalParams, 
    NumericVector & fitness,
    NumericVector & f,
    IntegerVector & trial,
    Function & fun,
    NumericMatrix & Foods,
    const NumericVector & lb,
    const NumericVector & ub
  ) {
  
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
    if (solution[param2change] < lb.at(param2change))
      solution[param2change] = lb.at(param2change);
    if (solution[param2change] > ub.at(param2change))
      solution[param2change] = ub.at(param2change);
    
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

void CalculateProbabilities(
    NumericMatrix & Foods,
    NumericVector & fitness,
    NumericVector & prob
  ) {
  double maxfit = fitness[0];
  for (int i=0;i<Foods.nrow();i++) 
    if (fitness[i] > maxfit) maxfit = fitness[i];
  
  for (int i=0;i<Foods.nrow();i++) 
    prob[i] = (0.9*((fitness[i] + 1e-40)/(maxfit + 1e-40))) + 0.1;
  
  return;
  
}

void SendOnlookerBees(
    double & GlobalMin,
    NumericVector & GlobalParams, 
    NumericVector & fitness,
    NumericVector & f,
    IntegerVector & trial,
    NumericVector & prob,
    Function & fun,
    NumericMatrix & Foods,
    const NumericVector & lb,
    const NumericVector & ub) {
  
  int param2change, neighbour;
  double ObjValSol, FitnessSol;
  NumericVector solution(Foods.ncol());
  
  int t = 0, i=0;
  while (t < Foods.nrow()) {
    // Randomly choose a food source
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
      if (solution[param2change] < lb.at(param2change))
        solution[param2change] = lb.at(param2change);
      
      if (solution[param2change] > ub.at(param2change))
        solution[param2change] = ub.at(param2change);
      
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
    NumericVector & fitness,
    NumericVector & f,
    IntegerVector & trial,
    NumericVector & prob,
    Function & fun,
    NumericMatrix & Foods,
    const NumericVector & lb,
    const NumericVector & ub,
    int limit
) {
  
  int maxtrialindex = 0;
  for(int i=0;i<Foods.nrow();i++) {
    if (trial[i] > trial[maxtrialindex])
      maxtrialindex = i;
  }
  
  // If it has reach the max, then init again
  if  (trial[maxtrialindex]>=limit)
    init(maxtrialindex, fitness, f, trial, fun, Foods, lb, ub);
  
  return;
}

// [[Rcpp::export]]
List abc_cpp_(
    NumericVector & par, 
    Function & fn,
    const NumericVector & lb,
    const NumericVector & ub,
    int FoodNumber   = 20, 
    int limit        = 100,
    int maxCycle     = 1000,
    int criter       = 50 // , double tol=1e-10
) {
  
  // Initialize:
  NumericMatrix Foods(FoodNumber, par.size());
  NumericVector f(FoodNumber), fitness(FoodNumber), prob(FoodNumber);
  IntegerVector trial(FoodNumber);
  
  NumericMatrix ans(maxCycle, par.size());

  // Initializing
  NumericVector GlobalParams = clone(par);
  double GlobalMin = as<double>(fn(GlobalParams));

  // Should be distributed equally
  for (int i=0;i<FoodNumber;i++) {
    
    for (int j=0;j<par.size();j++)
      Foods.at(i,j) = lb.at(j) + (ub.at(j) - lb.at(j))/(FoodNumber-1.0)*i;
    
    // Checking if it is defined
    f[i] = as<double>(fn(Foods(i,_)));
    if (NumericVector::is_na(f[i]))
      stop("Undefined value for -fn-. Check the function's domain.");
    
    fitness[i] = CalculateFitness(f[i]);
    trial[i]   = 0;
  }
  
  int unchanged = 0;
  MemorizeBetsSource(GlobalMin, GlobalParams, f, Foods, unchanged);
  ans(0,_) = GlobalParams;
  unchanged = 0;

  // double OldGlobalMin = GlobalMin;
  int i=0;
  while (++i<maxCycle) {
    
    SendEmployedBees(GlobalMin, GlobalParams, fitness, f, trial, fn, Foods, lb, ub);
    CalculateProbabilities(Foods, fitness, prob);
    SendOnlookerBees(GlobalMin, GlobalParams, fitness, f, trial, prob, fn, Foods, lb, ub);
    MemorizeBetsSource(GlobalMin, GlobalParams, f, Foods, unchanged);
    
    // Storing and breaking
    ans(i,_) = GlobalParams;
    if (unchanged>=criter) break;
    
    // If not, then send them out again
    SendScoutBees(fitness, f, trial, prob, fn, Foods, lb, ub, limit);
  }
  
  List obj = List::create(
    _["Foods"]  = Foods,
    _["f"]      = f,
    _["fitness"]= fitness,
    _["trial"]  = trial,
    _["value"]  = GlobalMin,
    _["par"]    = GlobalParams,
    _["counts"] = i,
    _["hist"]   = ans(Range(0,i),_)
  );
  
  obj.attr("class") = "abc_answer";
  
  return obj;
}


/***R
fun <- function(x) {
  -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2))
}

library(microbenchmark)
library(ABCoptim)
microbenchmark(
  abc_cpp(c(1,1),fun, ub = 5,lb = -5, criter = 20)[-8],
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