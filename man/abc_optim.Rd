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
  \item{par}{Initial values for the parameters to be optimized over}
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
  \item{criter}{Stop criteria (numer of unchanged results) until stopping}
}
\details{
  This is an implementation of Karaboga (2005) ABC optimization algorithm. It was developed upon the basic version programmed in \code{C} and distributed at the algorithm's official website (see the references).
  
  By default, lower and upper bounds are set as \code{+-Inf}. This last thing is just conceptual as all infinite bounds are replaced by \code{.Machine$double.xmax*1e-10} (still a pretty big number).
  
  If \code{D} (the number of parameters to be optimzed) is greater than one, then \code{lb} and \code{ub} can be either scalars (assuming that all the parameters share the same boundaries) or vectors (the parameters have different boundaries each other).
}
\value{
  A list containing the optimized parameters (\code{$par}), the value of the 
  function (\code{$value}) and the number of iterations taken to reach the optimum (\code{$counts}).
}
\author{George Vega Yon \email{g.vegayon@gmail.com}}
\references{
D. Karaboga, \emph{An Idea based on Honey Bee Swarm for Numerical Optimization}, tech. report TR06,Erciyes University, Engineering Faculty, Computer Engineering Department, 2005 \url{http://mf.erciyes.edu.tr/abc/pub/tr06_2005.pdf}

Artificial Bee Colony (ABC) Algorithm (website) \url{http://mf.erciyes.edu.tr/abc/index.htm}

Basic version of the algorithm implemented in \code{C} (ABC's official website) \url{http://mf.erciyes.edu.tr/abc/form.aspx}
}
\keyword{ optimization }

\examples{
# EXAMPLE 1: The minimum is at (pi,pi)
fun <- function(x) {
  -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2))
}

abc_optim(rep(0,2), fun, lb=-20, ub=20, criter=100)

# EXAMPLE 2: global minimum at about (-15.81515)
fw <- function (x)
  10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80

abc_optim(50, fw, lb=-100, ub=100, criter=100)

EXAMPLE 3: 5D sphere, global minimum at about (0,0,0,0,0)
fs <- function(x) sum(x^2)

abc_optim(rep(10,5), fs, lb=-100, ub=100, criter=200)
}