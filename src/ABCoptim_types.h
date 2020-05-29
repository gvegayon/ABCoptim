#include <Rcpp.h>

#ifndef H_ABCOPTIM 
#define H_ABCOPTIM 1
class TempParams {
public:
  unsigned int param2change = 0u;
  unsigned int neighbour = 0u;
  double ObjValSol = 0.0;
  double FitnessSol = 0.0;
  unsigned int t = 0;
  unsigned int i = 0;
  Rcpp::NumericVector solution;
  TempParams(unsigned int n) : solution(n) {};
  TempParams() : solution(0u) {};
  ~TempParams() {};
};

#endif