[![Downloads](http://cranlogs.r-pkg.org/badges/ABCoptim)](http://cran.rstudio.com/web/packages/ABCoptim/index.html) [![](http://cranlogs.r-pkg.org/badges/grand-total/ABCoptim)](http://cran.rstudio.com/web/packages/ABCoptim/index.html) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/gvegayon/ABCoptim?branch=master&svg=true)](https://ci.appveyor.com/project/gvegayon/ABCoptim) [![Travis-CI Build Status](https://travis-ci.org/gvegayon/ABCoptim.svg?branch=master)](https://travis-ci.org/gvegayon/ABCoptim) [![Coverage Status](https://img.shields.io/codecov/c/github/gvegayon/ABCoptim/master.svg)](https://codecov.io/github/gvegayon/ABCoptim?branch=master)

ABCoptim (beta)
===============

An implementation of the Artificial Bee Colony (ABC) Algorithm (R-package)

This is an implementation of Karaboga (2005) ABC optimization algorithm. It was developed upon the basic version programmed in *C* and distributed at the algorithm's official website (see the references).

Please consider that this version is in alpha state of development, thus any evident (precision) error should be blaimed to the package author (not to the algorithm itself).

Example
=======

``` r
library(ABCoptim)

# Function to optimize. Min at (pi,pi)
fun <- function(x) {
  -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2))
}

# Since it is stochastic, we need to set a seed to get the same
# results.
set.seed(123)

# Finding the minimum
abc_optim(rep(10,2), fun, lb=-20, ub=20, criter=200)
```

    ## An object of class -abc_answer- (Artificial Bee Colony Optim.):List of 8
    ##  $ Foods  : num [1:20, 1:2] 3.14 3.14 3.14 3.14 3.14 ...
    ##  $ f      :function (x)  
    ##   ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 4 8 6 1 8 1 4 6
    ##   .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x2574350> 
    ##  $ fitness: num [1:20] 2 2 2 2 2 ...
    ##  $ trial  : num [1:20] 15 29 19 19 16 10 14 39 19 4 ...
    ##  $ value  : num -1
    ##  $ par    : num [1:2] 3.14 3.14
    ##  $ counts : Named num 484
    ##   ..- attr(*, "names")= chr "function"
    ##  $ hist   : num [1:484, 1:2] 3.16 3.16 3.16 3.16 3.16 ...
    ##  - attr(*, "class")= chr "abc_answer"

Authors
=======

George G. Vega Yon \[aut\]

Enyelbert Muñoz \[cnt\]

References
==========

D. Karaboga, *An Idea based on Honey Bee Swarm for Numerical Optimization*, tech. report TR06,Erciyes University, Engineering Faculty, Computer Engineering Department, 2005 <http://mf.erciyes.edu.tr/abc/pub/tr06_2005.pdf>

Artificial Bee Colony (ABC) Algorithm (website) <http://mf.erciyes.edu.tr/abc/index.htm>

Basic version of the algorithm implemented in 'C' (ABC's official website) <http://mf.erciyes.edu.tr/abc/form.aspx>
