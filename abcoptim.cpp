#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double abcOptim_calculateFitness(double x) {  
	if (x >= 0) return 1/(x + 1);
	else return(1 + fabs(x));
}

// The best food source is memorized
// [[Rcpp::export]]
void abcOptim_memorizeBestSource(
	List foods,
	List global_min	
) {

	// Poiting to the foods position en valfun
	NumericMatrix foods_pos = foods["pos"];
	NumericVector foods_val = foods["val"];
	int nfoods = foods_val.size();
	int nparam = foods_pos.ncol();

	/* Pointing to the globmin pos and valfun */
	NumericMatrix global_min_pos = global_min["pos"];
	NumericVector global_min_val = global_min["val"];
	NumericVector global_min_per = global_min["per"];

	for(int i=0;i<nfoods;i++)
		if (foods_val[i] < global_min_val[0]) {
      
			/* Found a new global minimum */
			global_min_val[0] = foods_val[i];
			for(int j=0;j<nparam;j++)
				global_min_pos(0,j) = foods_pos(i,j);
        
		}
		else global_min_per[0]+=1;

	return ;
}

// [[Rcpp::export]]
NumericVector abcOptim_fun(Function targetfun, NumericVector par) {
	return targetfun(par);
}

// [[Rcpp::export]]
List abcOptim_gen(int nfoods,int nparam) {
	NumericMatrix pos(nfoods,nparam);
	NumericVector val(nfoods);
	NumericVector per(nfoods);

	return List::create(
		_["pos"] = pos,
		_["val"] = val,
		_["per"] = per
	);
} 

// [[Rcpp::export]]
void abcOptim_init(
	List Foods,
	bool optiinteger,
	NumericVector lb, NumericVector ub,
	Function targetfun
) {
	int i,j;

	/* Getting the food number */
	NumericMatrix foods_pos = Foods["pos"];
	NumericVector foods_val = Foods["val"];
	int nfoods = foods_pos.nrow();
	int nparam = foods_pos.ncol();

	if (optiinteger) {
		for(i=0;i<nfoods;i++)
			foods_pos(i,_) = runif(nparam) > .5;

			/* Getting the first value */
			foods_val[i] = abcOptim_fun(targetfun, foods_pos(i,_))[0];
	}
	else {
		for(i=0;i<nfoods;i++) {
			/* Setting the positions */
			for(j=0;j<nparam;j++)
				foods_pos(i,j) = runif(1, lb[j], ub[j])[0];
			
			/* Getting the first value */
			foods_val[i] = abcOptim_fun(targetfun, foods_pos(i,_))[0];
		}

	}

	return ;
}

// [[Rcpp::export]]
void abcOptim_sendEmployedBees(List Foods, bool optiinteger) {
	/* Poiting */
	NumericMatrix foods_pos = Foods["pos"];
	NumericVector foods_val = Foods["val"];

	NumericMatrix foods_pos0 = clone(foods_pos);
	NumericVector foods_val0 = clone(foods_val);

	int nparam = foods_pos.ncol();
	int nfoods = foods_val.size();

	NumericVector r = runif(nfoods);

	int param2change = floor(r[0]*nparam);

	/* Random solution */
	IntegerVector neighbour(nfoods);
	for(int i=0;i<nfoods;i++) 
		while(neighbour[i]==i) neighbour[i] = floor(r[i]*nfoods);

	if (optiinteger) foods_pos0[_,param2change] = r > .5;

	return ;
}

/*** R
foods <- list(
  pos=
    rbind(
      c(0,1),c(2,3),c(4,5)
      ),
  val=c(-1,.1,1)
)
foods
global_min <- list(pos=matrix(c(1,1),nrow=1),val=.15,per=0)
global_min
abcOptim_memorizeBestSource(foods,global_min)
global_min
*/

/* // [[Rcpp::export]] */



