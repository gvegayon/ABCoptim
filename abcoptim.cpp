#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double CalculateFitness(double x) {  
	if (x >= 0) return 1/(x + 1);
	else return(1 + fabs(x));
}

// The best food source is memorized
// [[Rcpp::export]]
List MemorizeBestSource(
	List foods,
	List global_min	
) {

	// Poiting to the foods position en valfun
	NumericMatrix foods_pos = foods["pos"];
	NumericVector foods_val = foods["val"];
	int nfoods = foods_val.size();
	int nparam = foods_pos.ncol();

	/* Pointing to the globmin pos and valfun */
	NumericVector global_min_pos = global_min["pos"];
	double global_min_val = global_min["val"];

	/* New potential value */
	NumericVector new_global_min_pos(global_min_pos);
	double new_global_min_val = global_min_val;
	int new_persistance = global_min["per"];

	for(int i=0;i<nfoods;i++)
		if (foods_val[i] < new_global_min_val) {
      
			/* Found a new global minimum */
			new_global_min_val = foods_val[i];
			for(int j=0;j<nparam;j++)
				new_global_min_pos[j] = foods_pos(i,j);
        
		}
		else new_persistance+=1;

	/* Returning output */
	return List::create(
		_["pos"]=new_global_min_pos,
		_["val"]=new_global_min_val,
		_["per"]=new_persistance
	);
}

// [[Rcpp::export]]
NumericVector fun(Function targetfun, NumericVector par) {
	return targetfun(par);
} 

// [[Rcpp::export]]
List init(
	List Foods,
	int index, bool firstinit, bool optiinteger,
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
			foods_val[i] = fun(targetfun, foods_pos(i,_))[0];
	}
	else {
		for(i=0;i<nfoods;i++) {
			/* Setting the positions */
			for(j=0;j<nparam;j++)
				foods_pos(i,j) = runif(1, lb[j], ub[j])[0];
			
			/* Getting the first value */
			foods_val[i] = fun(targetfun, foods_pos(i,_))[0];
		}

	}

	return List::create(
		_["pos"] = foods_pos,
		_["val"] = foods_val
	);
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
global_min <- list(pos=c(1,1),val=.15,per=0)
MemorizeBestSource(foods,global_min)
*/

/* // [[Rcpp::export]] */



