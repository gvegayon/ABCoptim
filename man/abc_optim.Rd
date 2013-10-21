\name{abc_optim}
\alias{abc_optim}
\title{Optimization through ABC algorithm}
\description{Optimizes through the ABC algorithm}
\usage{
abc_optim(par, fn, D=length(par), ..., NP=40, FoodNumber = NP/2,
    lb=-Inf, ub=+Inf, limit=100, maxCycle=1000, optiinteger=FALSE,
    criter=50)
}
\arguments{
  \item{par}{Initioal values for the parameters to be optimized over}
  \item{fn}{A function to be minimized, with first argument of the vector of parameters over which minimization is to take place. It should return a scalar result.}
  \item{D}{Number of parameters to be optimized.}
  \item{...}{Further arguments to be passed to 'fn'.}
  \item{NP}{Number of bees.}
  \item{FoodNumber}{Number of food sources to exploit.}
  \item{lb}{Lower bound of the parameters to be optimized.}
  \item{ub}{Upper bound of the parameters to be optimized.}
  \item{limit}{Limit of a food source.}
  \item{maxCycle}{Maximum number of iterations.}
  \item{optiinteger}{Whether to optimize binary parameters or not.}
  \item{criter}{Stop criteria (numer of unchanged results}
}
\value{
  A list containing the optimized parameters, the value of the 
  function and the number of iterations taken to reach the optimum.
}
\author{George Vega Yon \email{g.vegayon@gmail.com}}

\keyword{ optimization }
