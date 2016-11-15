# Changes in ABCoptim version 0.20 (2016-11-15)

* Added a `NEWS.md` file to track changes to the package.

* New function `abc_cpp` written with Rcpp (roughly 50%-100% faster than
  `abc_optim`).
  
* Both `abc_cpp` and `abc_optim` return objects of class `abc_answer`.

* Objects of class `abc_answer` return the trace of the global optimums.



