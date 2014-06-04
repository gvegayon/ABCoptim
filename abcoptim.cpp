#include <Rcpp.h>

using namespace Rcpp;

/* Calculates runif */
// [[Rcpp::export]]
double abc_unif(double minval=0.0, double maxval=1.0)
{
  return ((double)rand()/(double)RAND_MAX)*(maxval - minval) + minval;
}

// [[Rcpp::export]]
double abc_calc_fit(double x) {  
	if (x >= 0) return 1/(x + 1);
	else return(1 + fabs(x));
}

// The best food source is memorized
// [[Rcpp::export]]
void abc_mem_best_src(List Foods)
{
	// Poiting to the foods position en valfun
	NumericMatrix foods_pos = Foods["pos"];
	NumericVector foods_val = Foods["val"];
	int nfoods = foods_val.size();
	int nparam = foods_pos.ncol();

	/* Pointing to the globmin pos and valfun */
	NumericVector global_pos  = Foods["GlobalParams"];
	NumericVector global_val  = Foods["GlobalMin"];
	IntegerVector persistance = Foods["persistance"];

	for(int i=0;i<nfoods;i++)
		if (foods_val[i] < global_val[0]) {
      
			/* Found a new global minimum */
			global_val[0] = foods_val[i];
			for(int j=0;j<nparam;j++)
				global_pos[j] = foods_pos(i,j);
        
		}
		else persistance[0]+=1;

	return ;
}

// [[Rcpp::export]]
double abc_fun(Function targetfun, NumericVector par) {
	return as<double>(targetfun(par));
}

// [[Rcpp::export]]
List abc_initialize(
  int nfoods,
  int nparam,
  Function objfun,
  NumericVector lb,
  NumericVector ub,
  int limit
) 
{
  NumericVector fitness(nfoods);
	NumericMatrix pos(nfoods,nparam);
	NumericVector val(nfoods);
	NumericVector prob(nfoods);
  NumericVector trials(nfoods);
  NumericVector GlobPars(nparam);
  
  for(int i=0;i<nfoods;i++)
  {
    for(int j=0;j<nparam;j++)
      pos(i,j) = abc_unif(lb[j],ub[j]);
      
    val[i]     = as<double>(objfun(pos(i,_)));
    fitness[i] = abc_calc_fit(val[i]);
  }
  

	return List::create(
    _["fitness"]      = fitness,  /* Fitness of each food source */
    _["D"]            = nparam,   /* Number of parameters */
    _["FoodNumber"]   = nfoods,   /* Number of food sources */
		_["pos"]          = pos,      /* Current possition of each food source */
    _["prob"]         = prob,
		_["val"]          = val,      /* Current value of each food source */
    _["GlobalMin"]    = val[0],   /* Current Global Minimum */
    _["GlobalParams"] = pos(0,_), /* Current optimal params */
    _["trials"]       = trials,   /* Number of trials on the source */
    _["persistance"]  = 0,        /* N of iter without improvement */
    _["lb"]           = lb,       /* Lower bound */
    _["ub"]           = ub,       /* Upper bound */
    _["limit"]        = limit     /* Maximum number of attemps */
	);
} 

// [[Rcpp::export]]
void abc_init(
  int index,
	List Foods,
	Function objfun
) {
	int i,j;

	/* Getting the food number */
	NumericMatrix foods_pos = Foods["pos"];
	NumericVector foods_val = Foods["val"];
	int nparam = foods_pos.ncol();
  NumericVector lb = Foods["lb"];
  NumericVector ub = Foods["ub"];

	/* Setting the positions */
	for(j=0;j<nparam;j++)
		foods_pos(index,j) = abc_unif(lb[j],ub[j]);
	
	/* Getting the first value */
	foods_val[index] = abc_fun(objfun, foods_pos(index,_));

	return ;
}

// [[Rcpp::export]]
void abc_send_empl_bees(List Foods) {
	/* Poiting */
	NumericMatrix foods_pos = Foods["pos"];
	NumericVector foods_val = Foods["val"];

	NumericMatrix foods_pos0=clone(foods_pos);
	NumericVector foods_val0=clone(foods_val);

	int nparam = Foods["D"];
	int nfoods = Foods["FoodNumber"];

	NumericVector r = runif(nfoods+1);

	int param2change = floor(r[0]*nparam);

	/* Random solution */
	IntegerVector neighbour(nfoods);
	for(int i=0;i<nfoods;i++) 
		while(neighbour[i]==i) neighbour[i] = floor(r[i+1]*nfoods);

  /* Changin the values */
	for(int i=0;i<nfoods;i++) 
		foods_pos0(i,param2change) =
			foods_pos0(i,param2change) +
			(foods_pos0(i,param2change) - 
			foods_pos0(neighbour[i],param2change)*(r[i+1]-.5)*2);

	return ;
}

/* A food source is chosen with the probability which is proportioal to its 
quality. Different schemes can be used to calculate the probability values For
example prob(i)=fitness(i)/sum(fitness) or in a way used in the metot below
prob(i)=a*fitness(i)/max(fitness)+b probability values are calculated by using
fitness values and normalized by dividing maximum fitness value
*/
// [[Rcpp::export]]
double abc_calc_prob(List Foods) {
  NumericVector fitness=Foods["fitness"];
  int FoodNumber = Foods["FoodNumber"];
  double maxfit  = fitness[0];
  
  for(int i=1;i<FoodNumber;i++) 
    if (fitness[i] > maxfit) maxfit = fitness[i];
  
  double prob = 0.0;
  for(int i=0;i<FoodNumber;i++)
    prob += (0.9*(fitness[i]/maxfit) + .1);
    
  return prob;
}

/* Pointer Function Definition */
/*type double (*ptrFun)(NumericVector x);*/
void abc_send_onlooker_bees(List Foods, Function objfun)
{
  // Onlooker Bee phase
  int i = 0;
  int t = 0;
  
  int FoodNumber          = Foods["FoodNumber"];
  int D                   = Foods["D"];
  NumericVector foods_val = Foods["val"];
  NumericMatrix foods_pos = Foods["pos"];
  NumericVector foods_trl = Foods["trials"];
  NumericVector foods_prb = Foods["prob"];
  // Rprintf("ok\n");
  NumericVector foods_fit = Foods["fitness"];
  
  NumericVector lb = Foods["lb"];
  NumericVector ub = Foods["ub"];
  
  double r=0.0;
  while(t < FoodNumber)
  {
    r = abc_unif();
    // choose a food source depending on its probability to be chosen
    if (r < foods_prb[i])
    {
      t += 1;
      r = abc_unif();
      
      // The parameter to be changed is determined randomly
      int param2change = floor(r*D);
      
      /* A randomly chosen solution is used in producing a mutant solution of
      the solution i Randomly selected solution must be different from the
      solution i*/        
      int neighbour = i;
      while(neighbour==i)
        neighbour = floor(abc_unif()*FoodNumber);

      NumericVector solution = foods_pos(i,_);
      
      /* v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
      r = abc_unif();
      
      /* if (optiinteger) solution[param2change] <<- r > .5
      else 
      {*/
      solution[param2change] =
        solution[param2change]+
        (solution[param2change]-foods_pos[neighbour,param2change])*(r-0.5)*2;
      
      /* if generated parameter value is out of boundaries, it is shifted onto
      the boundaries*/
      if (solution[param2change]<lb[param2change]) 
        solution[param2change] = lb[param2change];
      
      if (solution[param2change]>ub[param2change]) 
        solution[param2change] = ub[param2change];
        
      /*}*/
      
      double ObjValSol  = as<double>(objfun(solution));
      double FitnessSol = abc_calc_fit(ObjValSol);
      
      /* a greedy selection is applied between the current solution i and its
      mutant*/
      if (FitnessSol>foods_fit[i])
      {
        /* If the mutant solution is better than the current solution i,
        replace the solution with the mutant and reset the trial counter of
        solution i*/
        foods_trl[i]   = 0;
        foods_pos(i,_) = clone(solution);
        foods_val[i]   = ObjValSol;
        foods_fit[i]   = FitnessSol;
      } /* if the solution i can not be improved, increase its trial counter*/
      else foods_trl[i] += 1;
    }
    if (i<FoodNumber) i+= 1;
    else i = 0;
    
    // end of onlooker bee phase
  }
}

/* determine the food sources whose trial counter exceeds the "limit" value.
In Basic ABC, only one scout is allowed to occur in each cycle*/
void abc_send_scout_bees(List Foods, Function objfun)
{
  NumericVector trials = Foods["trials"];
  int FoodNumber       = Foods["FoodNumber"];
  int maxtrialindex    = 1;
  int limit            = Foods["limit"];
  for (int i=0;i<FoodNumber;i++)
    if (trials[i] > trials[maxtrialindex]) maxtrialindex = i;
  
  if (trials[maxtrialindex] >= limit)
    abc_init(maxtrialindex, Foods, objfun);
  
  
  return;
}

// [[Rcpp::export]]
List abc_optimCpp(
  NumericVector par, Function objfun,
  NumericVector lb = NumericVector::create(1,-DBL_MAX),
  NumericVector ub = NumericVector::create(1,DBL_MAX),
  int FoodNumber   = 20,
  int limit        = 100,
  int maxCycle     = 1000,
  int criter       = 50
  )
{

  /* Number of parameters  */
  int nparam = par.size();
  
  /* Boundaries */
  NumericVector lb0(nparam,lb[0]);
  NumericVector ub0(nparam,ub[0]);
  if (lb.size() == nparam) lb0 = clone(lb);
  if (ub.size() == nparam) ub0 = clone(ub);
  
  List Foods = abc_initialize(FoodNumber, nparam, objfun, lb0, ub0, limit);
  
  /* Memorizes the initial sol */
  abc_mem_best_src(Foods);
  NumericVector persistance = Foods["persistance"];
  /* Start! */
  int iter = 0;
  while (++iter < maxCycle)
  {
    
    abc_send_empl_bees(Foods);
    abc_calc_prob(Foods);
    Rprintf("calc_prob... Iter %d\n",iter);
    abc_send_onlooker_bees(Foods, objfun);
    Rprintf("send_onlooker... Iter %d\n",iter);
    abc_mem_best_src(Foods);
    if (persistance[0] > criter) break;
    abc_send_scout_bees(Foods, objfun);
  }

  return(
    List::create(
      _["par"]=Foods["GlobalParams"],
      _["value"]=Foods["GlobalMinimum"],
      _["counts"]=iter
      )
    );
}

/*** R

fun <- function(x) {
  -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2))
}

abc_optimCpp(rep(0,2), fun, lb=-20, ub=20, criter=100)


*/

/* // [[Rcpp::export]] */



