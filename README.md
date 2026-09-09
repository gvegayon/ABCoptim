

[![R-CMD-check-final](https://github.com/gvegayon/ABCoptim/actions/workflows/r.yml/badge.svg)](https://github.com/gvegayon/ABCoptim/actions/workflows/r.yml)
[![pkgdown](https://github.com/gvegayon/ABCoptim/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/gvegayon/ABCoptim/actions/workflows/pkgdown.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/ABCoptim)](https://cran.r-project.org/package=ABCoptim)
[![Downloads](https://cranlogs.r-pkg.org/badges/ABCoptim)](http://cran.rstudio.com/package=ABCoptim)
[![](https://cranlogs.r-pkg.org/badges/grand-total/ABCoptim)](http://cran.rstudio.com/package=ABCoptim) 
[![codecov](https://codecov.io/gh/gvegayon/ABCoptim/branch/master/graph/badge.svg)](https://codecov.io/gh/gvegayon/ABCoptim)
[![DOI](https://zenodo.org/badge/13732591.svg)](https://zenodo.org/badge/latestdoi/13732591)
[![Sponsor](https://img.shields.io/badge/-Sponsor-fafbfc?logo=GitHub%20Sponsors)](https://github.com/sponsors/gvegayon)

# ABCoptim <img src="man/figures/logo.png" data-fig-align="right" height="128" />

<!-- how-to-cite -->

> \[!NOTE\] **How to cite ABCoptim.** If you use **ABCoptim** in
> published work, please cite it:
>
> Vega Yon G, Muñoz E. *ABCoptim: Implementation of Artificial Bee
> Colony (ABC) Optimization*.
> doi:[10.32614/CRAN.package.ABCoptim](https://doi.org/10.32614/CRAN.package.ABCoptim)
>
> Run `citation("ABCoptim")` in R for the BibTeX entry.
> <!-- how-to-cite -->

This is an implementation of Karaboga (2005) ABC optimization algorithm.
It was developed upon the basic version programmed in *C* and
distributed at the algorithm’s official website (see the references).

Any evident (precision) error should be blamed to the package author
(not to the algorithm itself).

# Example

``` r
library(ABCoptim)
```

    Using ABCoptim in your research? Please cite it: citation("ABCoptim")

``` r
# Function to optimize. Min at (pi,pi)
fun <- function(x) {
  -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2))
}

# Since it is stochastic, we need to set a seed to get the same
# results.
set.seed(123)

# Finding the minimum
ans <- abc_optim(rep(10,2), fun, lb=-20, ub=20, criter=200)
ans
```


     An object of class -abc_answer- (Artificial Bee Colony Optim.):
     par:
        x[1]:  3.141593
        x[2]:  3.141593

     value:
              -1.000000

     counts:
               457

``` r
plot(ans)
```

![](man/figures/example1-1.png)

# References

D. Karaboga, *An Idea based on Honey Bee Swarm for Numerical
Optimization*, tech. report TR06,Erciyes University, Engineering
Faculty, Computer Engineering Department, 2005
http://mf.erciyes.edu.tr/abc/pub/tr06_2005.pdf

Artificial Bee Colony (ABC) Algorithm (website)
http://mf.erciyes.edu.tr/abc/index.htm

Basic version of the algorithm implemented in ‘C’ (ABC’s official
website) http://mf.erciyes.edu.tr/abc/form.aspx

## Code of Conduct

Please note that the ABCoptim project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
