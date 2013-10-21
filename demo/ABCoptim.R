################################################################################
# Some examples of ABC algorithm
# Author: George G. Vega
################################################################################
pause <- function() {  
  invisible(readline("\nPress <return> to continue: ")) 
}
pause()

## Rosenbrock Banana function, global minimum at about (1,1)
# No so good example (for now)
fr <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
pause()

abc_optim(c(50,50), fr, lb=-100, ub=100, criter=500)


## "wild" function , global minimum at about -15.81515
# Very good example
pause()
fw <- function (x)
  10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80

abc_optim(50, fw, lb=-100, ub=100, criter=100)

## Griewank function, global minimum at 0
# Very good example
pause()

fg <- function(x) 
  sum(x*x)/4000-prod(cos(x/sqrt(1:2)))+1

abc_optim(50, fg, lb=-100, ub=100, criter=100)

# Rastrigin function, global minimum at (0,0)
# Another very good example
pause()
fra <- function(x)
  20 + x[1]^2 + x[2]^2 - 10*(cos(2*pi*x[1]) + cos(2*pi*x[2]))

abc_optim(rep(50,2), fra, lb=-100, ub=100, criter=100)
