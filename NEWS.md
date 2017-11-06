# Changes in ABCoptim version 0.15.0 (2017-11-05)

* `abc_cpp` now checks the domain of the function.

* New arguments `parscale` and `fnscale` added (see stats::optim)

* New `print` and `plot` method for `abc_answer` class object.


# Changes in ABCoptim version 0.14.0 (2016-11-16)

* Added a `NEWS.md` file to track changes to the package.

* New function `abc_cpp` written with Rcpp (roughly 50%-100% faster than
  `abc_optim`).
  
* Both `abc_cpp` and `abc_optim` return objects of class `abc_answer`.

* Objects of class `abc_answer` return the trace of the global optimums.


