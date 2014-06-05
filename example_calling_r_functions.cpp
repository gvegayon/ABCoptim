#include <Rinternals.h>
#include <Rcpp.h>

// [[Rcpp::export]]
double abc_call_fun(SEXP par, SEXP fn)
{
  SEXP R_fcall, ans;

  if(!isFunction(fn)) error("'fn' must be a function");
  R_fcall = PROTECT(lang2(fn, R_NilValue));
  
  SETCADR(R_fcall, par);
  ans=eval(R_fcall, R_GlobalEnv);
  UNPROTECT(1);
  return REAL(ans)[0];
}

using namespace Rcpp;

// [[Rcpp::export]]
double rcpp_call_fun(NumericVector par, Function fn)
{
  return as<double>(fn(par));
}



/*** R
fun <- function(x) {
  -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2))
}

x <- runif(2)
library(microbenchmark)
microbenchmark(
  abc_call_fun(x, fun), rcpp_call_fun(x,fun), times=10000
)


*/


