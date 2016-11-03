# Artificial bee colony

rm(list=ls())

# Control Parameters of ABC algorithm
NP          <- 40   # The number of colony size (employed bees+onlooker bees)
FoodNumber  <- NP/2 # The number of food sources equals the half of the colony size
limit       <- 50  # A food source which could not be improved through "limit" trials is abandoned by its employed bee
maxCycle    <- 1000  # The number of cycles for foraging {a stopping criteria}

# Problem specific variables
D           <- 4     # The number of parameters of the problem to be optimized
lb          <- 0 # Lower bound of the parameters.
ub          <- 1  # Upper bound of the parameters.
optiinteger <- TRUE

# How many times to see robustness
runtime     <- 1    # Algorithm can be run many times in order to see its robustness
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
GlobalMin   <- double(1)
GlobalParams<- double(D)
GlobalMins  <- double(runtime)
r           <- integer(1)

# Function
fun1 <- function(x) x[1]^2 - 2*x[2]

# Definicion de la funcion
fun3 <- function(x) {
  -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2))
}

fun <- function(lambda, x0 = 2, X = c(3,2,3,1), M = 1)
{
  norm(x0 - X%*%lambda) + exp(abs(sum(lambda > 0) - M))
}

#fun <- function(...) fun2(...)*-1


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
  for(i in seq(1,FoodNumber)) 
  {
    if (f[i] < GlobalMin) 
    {
      GlobalMin <<- f[i]
      for (j in seq(1,D))
      {
        GlobalParams[j]<<-Foods[i,j]
      }
    }
  }
}
# MemorizeBestSource()

# Variables are initialized in the range [lb,ub]. If each parameter has different range, use arrays lb[j], ub[j] instead of lb and ub 
# Counters of food sources are also initialized in this function

init <- function(index)
{
  for (j in 1:D)
  {
    r <- runif(1)
    if (optiinteger) Foods[index,j] <<- r > .5
    else Foods[index,j] <<- r*(ub-lb) + lb
    solution[j] <<- Foods[index,j]
  }
  f[index] <<- fun(solution)
  fitness[index] <<- CalculateFitness(f[index]);
  trial[index] <<- 0;
  
}
# init(2)

# All food sources are initialized
initial <- function() {
  for (i in 1:FoodNumber)
  {
    init(i)
  }
  GlobalMin <<- f[1]
  for (i in 1:D)
  {
    GlobalParams[i] <<- Foods[1,i]
  }
}

# initial()


SendEmployedBees <- function()
{
  for (i in 1:FoodNumber) 
  {
    # The parameter to be changed is determined randomly
    r <- runif(1)
    param2change <- floor(r*D)
    if (!param2change) param2change <- 1
    
    # A randomly chosen solution is used in producing a mutant solution of the solution i
    r <- runif(1)
    neighbour <- floor(r*FoodNumber)
    if (!neighbour) neighbour <- 1
    
    # Randomly selected solution must be different from the solution i
    while(neighbour==i)
    {
      r <- runif(1)
      neighbour <- floor(r*FoodNumber)
      if (!neighbour) neighbour <- 1
    }
    
    solution <<- Foods[i,]
    
    # v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) 
    r <- runif(1)
    
    if (optiinteger) solution[param2change] <<- r > 0.5
    else
    {
      solution[param2change] <<- Foods[i,param2change]+(Foods[i,param2change]-Foods[neighbour,param2change])*(r-0.5)*2
      # if generated parameter value is out of boundaries, it is shifted onto the boundaries
      if (solution[param2change]<lb) solution[param2change]<<-lb
      if (solution[param2change]>ub) solution[param2change]<<-ub
    }
    
    ObjValSol <<- fun(solution)
    FitnessSol <<- CalculateFitness(ObjValSol)
    
    # a greedy selection is applied between the current solution i and its mutant*/
    if (FitnessSol>fitness[i])
    {
      # If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
      trial[i] <<- 0;
      for(j in 1:D) Foods[i,j] <<- solution[j]
      f[i]<<- ObjValSol
      fitness[i]<<-FitnessSol
    }
    else
    { # the solution i can not be improved, increase its trial counter*/
          trial[i] <<- trial[i]+1
    }
    
  }
}


# A food source is chosen with the probability which is proportioal to its quality*/
# Different schemes can be used to calculate the probability values*/
# For example prob(i)=fitness(i)/sum(fitness)*/
# or in a way used in the metot below prob(i)=a*fitness(i)/max(fitness)+b*/
# probability values are calculated by using fitness values and normalized by dividing maximum fitness value*/
CalculateProbabilities <- function()
{
  maxfit <<- fitness[1]
  for (i in 1:FoodNumber) 
  {
    if (fitness[i] > maxfit) maxfit <<- fitness[i]
  }
  
  for (i in 1:FoodNumber)
  {
    prob[i] <<- (.9*(fitness[i]/maxfit)) + .1
  }
}

SendOnlookerBees <- function()
{
  # Onlooker Bee phase
  i <- 1
  t <- 0
  while (t < FoodNumber)
  {
    r = runif(1)
    if (r < prob[i]) # choose a food source depending on its probability to be chosen
    {
      t <- t + 1
      
      # The parameter to be changed is determined randomly
      r <- runif(1)
      param2change <- floor(r*D)
      if (!param2change) param2change <- 1
      
      # A randomly chosen solution is used in producing a mutant solution of the solution i
      r = runif(1)
      neighbour= floor(r*FoodNumber)
      if (!neighbour) neighbour <- 1
      neighbour <- (1:FoodNumber)[order(runif(FoodNumber))][1]
      
      #Randomly selected solution must be different from the solution i*/        
      while(neighbour==i)
      {
        r = runif(1)
        neighbour= floor(r*FoodNumber)
        if (!neighbour) neighbour <- 1
      }
      
      solution <<- Foods[i,]
      
      # v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
      r = runif(1)
      
      if (optiinteger) solution[param2change] <<- r > .5
      else 
      {
        # if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        if (solution[param2change]<lb) solution[param2change] <<- lb
        if (solution[param2change]>ub) solution[param2change] <<- ub
        solution[param2change] <<- Foods[i,param2change]+(Foods[i,param2change]-Foods[neighbour,param2change])*(r-0.5)*2
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

# determine the food sources whose trial counter exceeds the "limit" value. In Basic ABC, only one scout is allowed to occur in each cycle*/

SendScoutBees <- function() 
{
  maxtrialindex <- 1
  for (i in 1:FoodNumber)
  {
    if (trial[i] > trial[maxtrialindex]) maxtrialindex <- i
  }
  
  if (trial[maxtrialindex] >= limit) init(maxtrialindex)
}


# Debug de las funcionts
# debug(SendEmployedBees)
# debug(SendOnlookerBees)

# Main program of the ABC algorithm
for(run in 1:runtime)
{
  promedio <- 0
  
  initial()
  MemorizeBestSource()
  for (iter in 0:maxCycle)
  {
    SendEmployedBees() # ; print("enviados")
    CalculateProbabilities() # ; print("probab")
    SendOnlookerBees() # ; print("lookon")
    MemorizeBestSource() # ; print("mem")
    SendScoutBees() # ; print("scouts")
  }
  for(j in 1:D)
  {
    message(sprintf("GlobalParam[%d]: %f",j,GlobalParams[j]))
  }
  message(sprintf("%d. run: %f",run,GlobalMin))
  GlobalMins[run] <- GlobalMin
  promedio <- promedio + GlobalMin
}
