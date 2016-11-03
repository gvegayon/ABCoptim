#include <Rinternals.h>
#include <Rcpp.h>

// [[Rcpp::export]]
SEXP cppFuncall(SEXP par, SEXP fn)
{
  SEXP R_fcall, ans;

  if(!isFunction(fn)) error("'fn' must be a function");
  R_fcall = PROTECT(lang2(fn, R_NilValue));
  
  SETCADR(R_fcall, par);
  ans=eval(R_fcall, R_GlobalEnv);
  UNPROTECT(1);
  
  return ans;
}
 
using namespace Rcpp;

// [[Rcpp::export]]
SEXP RcppFuncall(NumericVector par, Function fn)
{
  return fn(par);
}



/*** R
# R function to be called
fun <- function(x) {
  w<-sapply(1:1e3,function(x) -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2)))
  1
}

# Input data
set.seed(3331)
x <- runif(1e3)

# Benchmarking
library(microbenchmark)
print(microbenchmark(
  cppFuncall(x, fun), RcppFuncall(x,fun), fun(x), times=1e3,
  unit="relative", control = list(warmup=100)
), signif=2)

# cppFuncall(letters, fun)
*/


