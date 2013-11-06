\name{ABCoptim-package}
\alias{ABCoptim-package}
\alias{ABCoptim}
\alias{abc}
\docType{package}
\title{
  An implementation of the Artificial Bee Colony (ABC) Algorithm
}
\description{
  This is an implementation of Karaboga (2005) ABC optimization algorithm. It was developed upon the basic version programmed in \code{C} and distributed at the algorithm's official website (see the references).
  
  Please consider that this version is in alpha state of development, thus any evident (precision) error should be blaimed to the package author (not to the algorithm itself)

  Please visit the project home for more information:
    \url{https://github.com/gvegayon/ABCoptim}.
} 
\details{
  \tabular{ll}{
    Package: \tab ABCoptim\cr
    Type: \tab Package\cr
    Version: \tab 0.13.11\cr
    Date: \tab 2013-11-05\cr
    License: \tab GPL version 2 or later\cr
  }
}
\author{George Vega Yon \email{g.vegayon@gmail.com} [aut],

Enyelbert Muñoz \email{enyeldoc2011@gmail.com} [cnt)}

\references{
D. Karaboga, \emph{An Idea based on Honey Bee Swarm for Numerical Optimization}, tech. report TR06,Erciyes University, Engineering Faculty, Computer Engineering Department, 2005 \url{http://mf.erciyes.edu.tr/abc/pub/tr06_2005.pdf}

Artificial Bee Colony (ABC) Algorithm (website) \url{http://mf.erciyes.edu.tr/abc/index.htm}

Basic version of the algorithm implemented in \code{C} (ABC's official website) \url{http://mf.erciyes.edu.tr/abc/form.aspx}
}
\examples{
  \dontrun{
    demo(ABCoptim) # Some functions...
  }
}
\keyword{ package }
