#' Artificial Bee Colony Optimization
#' 
#' Implements Karaboga (2005) Artificial Bee Colony (ABC) Optimization algorithm.
#' 
#' @param par Initial values for the parameters to be optimized over
#' @param fn A function to be minimized, with first argument of the vector of
#' parameters over which minimization is to take place. It should return a
#' scalar result.
#' @param ... Further arguments to be passed to 'fn'.
#' @param FoodNumber Number of food sources to exploit. Notice that the param
#' \code{NP} has been deprecated.
#' @param lb Lower bound of the parameters to be optimized.
#' @param ub Upper bound of the parameters to be optimized.
#' @param limit Limit of a food source.
#' @param maxCycle Maximum number of iterations.
#' @param optiinteger Whether to optimize binary parameters or not.
#' @param criter Stop criteria (numer of unchanged results) until stopping
#' @param parscale Numeric vector of length \code{length(par)}. Scale applied
#' to the parameters (see \code{\link[stats:optim]{optim}}).
#' @param fnscale Numeric scalar. Scale applied function. If \code{fnscale < 0},
#' then the problem becomes a maximization problem (see \code{\link[stats:optim]{optim}}).
#'
#' @details 
#' 
#' This implementation of the ABC algorithm was developed based on the basic
#' version written in \code{C} and published at the algorithm's official
#' website (see references).
#' 
#' \code{abc_optim} and \code{abc_cpp} are two different implementations of the
#' algorithm, the former using pure \code{R} code, and the later using \code{C++},
#' via the \pkg{Rcpp} package. Besides of the output, another important
#' difference between the two implementations is speed, with \code{abc_cpp}
#' showing between 50\% and 100\% faster performance.
#' 
#' Upper and Lower bounds (\code{ub}, \code{lb}) equal to infinite will be replaced
#' by either \code{.Machine$double.xmax} or \code{-.Machine$double.xmax}.
#' 
#' If \code{D} (the number of parameters to be optimzed) is greater than one,
#' then \code{lb} and \code{ub} can be either scalars (assuming that all the
#' parameters share the same boundaries) or vectors (the parameters have
#' different boundaries each other).
#' 
#' @return An list of class \code{abc_answer}, holding the following elements:
#' \item{Foods}{Numeric matrix. Last position of the bees.}
#' \item{f}{Numeric vector. Value of the function evaluated at each set of \code{Foods}.}
#' \item{fitness}{Numeric vector. Fitness of each \code{Foods}.}
#' \item{trial}{Integer vector. Number of trials at each \code{Foods}.}
#' \item{value}{Numeric scalar. Value of the function evaluated at the optimum.}
#' \item{par}{Numeric vector. Optimum found.}
#' \item{counts}{Integer scalar. Number of cycles.}
#' \item{hist}{Numeric matrix. Trace of the global optimums.}
#' 
#' @author George Vega Yon \email{g.vegayon@@gmail.com}
#' @references D. Karaboga, \emph{An Idea based on Honey Bee Swarm for
#' Numerical Optimization}, tech. report TR06,Erciyes University, Engineering
#' Faculty, Computer Engineering Department, 2005
#' \url{http://mf.erciyes.edu.tr/abc/pub/tr06_2005.pdf}
#' 
#' Artificial Bee Colony (ABC) Algorithm (website)
#' \url{http://mf.erciyes.edu.tr/abc/index.htm}
#' 
#' Basic version of the algorithm implemented in \code{C} (ABC's official
#' website) \url{http://mf.erciyes.edu.tr/abc/form.aspx}
#' @keywords optimization
#' @examples
#' 
#' # EXAMPLE 1: The minimum is at (pi,pi) ----------------------------------------
#' fun <- function(x) {
#'   -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2))
#' }
#' 
#' ans0 <- abc_optim(rep(0,2), fun, lb=-10, ub=10, criter=50)
#' ans0[c("par", "counts", "value")]
#' 
#' ans1 <- abc_cpp(rep(0,2), fun, lb=-10, ub=10, criter=50)
#' ans1[c("par", "counts", "value")]
#' 
#' # EXAMPLE 2: global minimum at about (-15.81515)
#' fw <- function (x)
#'   10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80
#' 
#' ans <- abc_optim(50, fw, lb=-100, ub=100, criter=100)
#' ans[c("par", "counts", "value")]
#' 
#' # EXAMPLE 3: 5D sphere, global minimum at about (0,0,0,0,0)
#' fs <- function(x) sum(x^2)
#' 
#' ans <- abc_optim(rep(10,5), fs, lb=-100, ub=100, criter=200)
#' ans[c("par", "counts", "value")]
#' 
#' @export abc_optim
#' @aliases abc_answer
abc_optim <- function(
  par,               # Vector de parametros a opti 
  fn,                # Funcion objetivo
  ...,               # Argumentos de la funcion (M, x0, X, etc.)
  FoodNumber  = 20,   # Fuentes de alimento 
  lb          = rep(-Inf, length(par)),        # Limite inferior de recorrido
  ub          = rep(+Inf, length(par)),        # Limite superior de recorrido
  limit       = 100,       # Limite con que se agota una fuente de alimento
  maxCycle    = 1000,   # Numero maximo de iteraciones 
  optiinteger = FALSE, # TRUE si es que queremos optimizar en [0,1] (binario)
  criter      = 50,
  parscale    = rep(1, length(par)),
  fnscale     = 1
)
{
  D <- length(par)
  
  # Checking limits
  if (length(lb) == 1 && length(par) > 1) lb <- rep(lb, D)
  if (length(ub) == 1 && length(par) > 1) ub <- rep(ub, D)

  lb[is.infinite(lb)] <- -.Machine$double.xmax*1e-10
  ub[is.infinite(ub)] <- .Machine$double.xmax*1e-10
  
  # Initial params
  Foods       <- matrix(double(FoodNumber*D), nrow=FoodNumber)
  f           <- double(FoodNumber)
  fitness     <- double(FoodNumber)
  trial       <- double(FoodNumber)
  prob        <- double(FoodNumber)
  solution    <- double(D)
  ObjValSol   <- double(1)
  FitnessSol  <- double(1)
  neighbour   <- integer(1)
  param2change<- integer(1)
  GlobalMin   <- fn(par, ...) # double(1)
  GlobalParams<- par #double(D)
  #GlobalMins  <- double(runtime)
  r           <- integer(1)

  # Fun
  fun <- function(par) fn(par/parscale, ...)/fnscale
  
  # Fitness function
  CalculateFitness <- function(fun)
  {
    if (fun >= 0) return(1/(fun + 1))
    else return(1 + abs(fun))
  }
  # CalculateFitness(f[1])
  
  # The best food source is memorized
  MemorizeBestSource <- function() 
  {
    oldGlobalMin <- GlobalMin
    for(i in seq(1,FoodNumber)) {
      if (f[i] < GlobalMin) {
        GlobalMin <<- f[i]
        
        # Replacing new group of parameters
        GlobalParams <<- Foods[i,]
      }
    }
    
    # Increasing persistance
    if (oldGlobalMin == GlobalMin) persistance <<- persistance + 1
    else persistance <<- 0
  }
  
  # Variables are initialized in the range [lb,ub]. If each parameter has
  # different range, use arrays lb[j], ub[j] instead of lb and ub 
  # Counters of food sources are also initialized in this function
  
  init <- function(index, ...) {
    if (optiinteger) Foods[index,] <<- runif(D) > .5
    else {
      Foods[index,] <<- sapply(1:D, function(k) runif(1,lb[k],ub[k]) )
    }

    solution <<- Foods[index,]
    
    f[index] <<- fun(solution)

    fitness[index] <<- CalculateFitness(f[index])
    trial[index] <<- 0
    
  }
  # init(2)
  
  # All food sources are initialized
  initial <- function() {
    # For the first initialization we set the bees at
    # specific places equaly distributed through the
    # bounds.
    Foods <<- 
      sapply(1:D, function(k) {
        seq(lb[k],ub[k],length.out=FoodNumber)
      }
      )
    
    for (i in 1:FoodNumber) {
      solution <<- Foods[i,]
      
      f[i] <<- fun(solution)
      
      fitness[i] <<- CalculateFitness(f[i])
      trial[i] <<- 0
    }
  }
  
  # initial()
  
  
  SendEmployedBees <- function() {
    for (i in 1:FoodNumber) {
      # The parameter to be changed is determined randomly
      param2change <- sample(1:D, 1) # floor(runif(1)*D) + 1 
      
      # A randomly chosen solution is used in producing a mutant solution of the solution i
      # Randomly selected solution must be different from the solution i
      neighbour <- i
      while(neighbour==i)
        neighbour <- sample(1:FoodNumber, 1) # floor(runif(1)*FoodNumber) + 1
      
      solution <<- Foods[i,]
      
      # v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) 

      if (optiinteger) solution[param2change] <<- runif(1) > 0.5
      else {
        solution[param2change] <<- 
          Foods[i,param2change]+
          (Foods[i,param2change]-Foods[neighbour,param2change])*(runif(1)-0.5)*2

        # if generated parameter value is out of boundaries, it is shifted onto the boundaries
        if (solution[param2change]<lb[param2change])
          solution[param2change]<<-lb[param2change]
        
        if (solution[param2change]>ub[param2change])
          solution[param2change]<<-ub[param2change]
      }
      
      ObjValSol <<- fun(solution)
      FitnessSol <<- CalculateFitness(ObjValSol)
      
      # a greedy selection is applied between the current solution i and its mutant*/
      if (FitnessSol>fitness[i]) {
        # If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
        trial[i] <<- 0;
        #for(j in 1:D) Foods[i,j] <<- solution[j]
        Foods[i,] <<- solution
        f[i]<<- ObjValSol
        fitness[i]<<-FitnessSol
      }
      else {
        # the solution i can not be improved, increase its trial counter*/
        trial[i] <<- trial[i]+1
      }
    }
  }
  
  
  # A food source is chosen with the probability which is proportioal to its quality*/
  # Different schemes can be used to calculate the probability values*/
  # For example prob(i)=fitness(i)/sum(fitness)*/
  # or in a way used in the metot below prob(i)=a*fitness(i)/max(fitness)+b*/
  # probability values are calculated by using fitness values and normalized by dividing maximum fitness value*/
  CalculateProbabilities <- function() {
    maxfit <- fitness[1]
    for (i in 1:FoodNumber) 
      if (fitness[i] > maxfit) maxfit <- fitness[i]
    
    prob <<- .9*(fitness/(maxfit+1e-20)) + .1
#     prob[is.nan(prob)]  <<- .1
  }
  
  SendOnlookerBees <- function()
  {
    # Onlooker Bee phase
    i <- 1
    t <- 0
    while (t < FoodNumber)
    {

      # choose a food source depending on its probability to be chosen
      if (runif(1) < prob[i]) {
        t <- t + 1

        # The parameter to be changed is determined randomly
        param2change <- sample(1:D, 1) # floor(runif(1)*D) + 1 
        
        # A randomly chosen solution is used in producing a mutant solution of the solution i
        #Randomly selected solution must be different from the solution i*/        
        neighbour <- i
        while(neighbour==i)
          neighbour <- sample(1:FoodNumber, 1) # floor(runif(1)*FoodNumber) + 1

        solution <<- Foods[i,]
        
        # v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */

        if (optiinteger) solution[param2change] <<- runif(1) > .5
        else 
        {
          solution[param2change] <<- 
            Foods[i,param2change]+
            (Foods[i,param2change]-Foods[neighbour,param2change])*(runif(1)-0.5)*2
          
          # if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
          if (solution[param2change]<lb[param2change]) 
            solution[param2change] <<- lb[param2change]
          
          if (solution[param2change]>ub[param2change]) 
            solution[param2change] <<- ub[param2change]
          
        }
        
        ObjValSol <<- fun(solution)
        FitnessSol <<- CalculateFitness(ObjValSol)
        
        # a greedy selection is applied between the current solution i and its mutant*/
        if (FitnessSol>fitness[i])
        {
          # If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
          trial[i] <<- 0
          Foods[i,] <<- solution
          
          f[i]<<-ObjValSol
          fitness[i]<<-FitnessSol
        } #if the solution i can not be improved, increase its trial counter*/
        else trial[i] <<- trial[i]+1
      }
      i <- i + 1
      if (i==FoodNumber) i <- 1
      # end of onlooker bee phase
    }
  }
  
  # determine the food sources whose trial counter exceeds the "limit" value.
  # In Basic ABC, only one scout is allowed to occur in each cycle*/
  
  SendScoutBees <- function() {
    maxtrialindex <- 1
    for (i in 1:FoodNumber) {
      if (trial[i] > trial[maxtrialindex]) maxtrialindex <- i
    }
    
    if (trial[maxtrialindex] >= limit) init(maxtrialindex)
  }
  
  persistance <- 0
  
  # Inicializa funcion
  initial()
  
  # Memoriza la primera mejor solucion
  MemorizeBestSource() 
  
  ans  <- matrix(0, ncol = D, nrow=maxCycle)
  iter <- 0
  # Comienza a iterar
  while ((iter <- iter + 1) < maxCycle)
  {
    SendEmployedBees()
    CalculateProbabilities()
    SendOnlookerBees() 
    MemorizeBestSource()
    
    # Storing parameter and breaking out
    ans[iter,] <- GlobalParams
    if (persistance > criter) break
    
    SendScoutBees()
  }

  return(
    structure(list(
      Foods   = Foods,
      f       = fn,
      fitness = fitness,
      trial   = trial,
      value   = fun(GlobalParams),
      par     = GlobalParams,
      counts  = c("function"=iter),
      hist    = ans[1:iter,,drop=FALSE]
      ), class="abc_answer"
    ))
  
}

#' @export
#' @param x An object of class \code{abc_answer}.
#' @rdname abc_optim
print.abc_answer <- function(x, ...) {
  cat("An object of class -abc_answer- (Artificial Bee Colony Optim.):")
  str(x)
  invisible(x)
}

# ################################################################################
# # Ejemplos
# ################################################################################
# 
# X <- c(3,2,3,1)
# 
# # Funcion de matching
# fun <- function(lambda, x0, X, M)
# {
#   norm((x0 - X)*lambda, type="2") + exp(abs(sum(lambda > 0) - M))
# }
# 
# # Mejor vecino para
# #  x0 = 2
# #  X  = c(3,2,3,1)
# #  M  = 1
# # El mejor resultado debe ser [0,1,0,0]
# x1 <- abc_optim(rep(0,4), fun, x0=2, X=X, M=1, lb=0, ub=1, optiinteger=T)
# x1
# 
# # Mejores dos vecinos para
# #  x0 = 3
# #  X  = c(3,2,3,1)
# #  M  = 2
# # El mejor resultado debe ser [1,0,1,0]
# x2 <- abc_optim(rep(0,4), fun, x0=3, X=X, M=2, lb=0, ub=1, optiinteger=T)
# x2
# 
# ################################################################################
# # Definicion de la funcion
# fun <- function(x) {
#   -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2))
# }
# 
# abc_optim(rep(0,2), fun, lb=-5, ub=5, criter=50)
# 
# optim(rep(0,2), fn=fun) #lower=-5,upper=5)
# 
# ################################################################################
# # Definicion de la funcion
# 
# fun <- function(x) {
#   -4+(x[1]^2 + x[2]^2)
# }
# 
# abc_optim(c(1,1), fn=fun, lb=-100000, ub=100000,criter=100)
# 
# ################################################################################
# # Definicion de la funcion
# 
# fun <- function(x) {
#   -(x^4 - 2*x^2 - 8)
# }
# 
# abc_optim(0, fn=fun, lb=-2, ub=2,criter=100)
# # 


# library(microbenchmark)
# const <- 2
# fun <- function(x) {
#   -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^const + (x[2] - pi)^const))
# }
# 
# microbenchmark(
#   ABC_R = abc_optim(rep(0,2), fun, lb=-20, ub=20, criter=20, maxCycle = 20),
#   ABC_CPP = abc_cpp(rep(0,2), fun, lb=-20, ub=20, criter=20, maxCycle = 20),
#   times=100
# )


#' @export
#' @rdname abc_optim
abc_cpp <- function(
  par,
  fn,
  ...,
  FoodNumber = 20,   # Fuentes de alimento 
  lb         = rep(-Inf, length(par)),        # Limite inferior de recorrido
  ub         = rep(+Inf, length(par)),        # Limite superior de recorrido
  limit      = 100,       # Limite con que se agota una fuente de alimento
  maxCycle   = 1000,   # Numero maximo de iteraciones 
  criter     = 50,
  parscale   = rep(1, length(par)),
  fnscale    = 1
) {
  
  # Checking limits
  if (length(lb)>0) lb <- rep(lb, length(par))
  if (length(ub)>0) ub <- rep(ub, length(par))
  
  lb[is.infinite(lb)] <- -(.Machine$double.xmax*1e-10)
  ub[is.infinite(ub)] <- +(.Machine$double.xmax*1e-10)
  
  fun <- function(par) fn(par/parscale, ...)/fnscale
  
  abc_cpp_(par, fun, lb, ub, FoodNumber, limit, maxCycle, criter)
  
}