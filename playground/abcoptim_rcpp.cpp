#include <ABCoptim.h>

using namespace Rcpp;

// [[Rcpp::export]]
double funRcpp(NumericVector x)
{
  return -cos(x[0])*cos(x[1])*exp(-(
    pow(x[0]-M_PI,2.0) + pow(x[1]-M_PI,2.0) ));
}

// [[Rcpp::export]]
List abc_optim2(
  NumericVector par,
  NumericVector lb,
  NumericVector ub,
  int FoodNumber   = 20,
  int limit        = 100,
  int maxCycle     = 1000,
  int criter       = 50)
{
 return(
   abc_optim_Cpp(
     par,&funRcpp, lb,ub,FoodNumber,limit,maxCycle,criter
     )
  );
}


/*** R

fun <- function(x) {
  -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2))
}

library(microbenchmark)
library(ABCoptim)

message("Funcion en C")
x1<-abc_optim2(par=rep(0,2), lb=-20, ub=20, criter=200)
x1

#microbenchmark(
#abc_optimCpp(par=rep(0,2), lb=-20, ub=20, criter=200),times=1000
#)

system.time(abc_optim2(par=rep(0,2), lb=-20, ub=20, criter=200))


***/ 





