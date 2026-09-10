- Use the tinyverse (not tidyverse) principle: keep dependencies to a minimum.
- Always document code using `roxygen2`.
- Whe adding new features/fixing bugs/improving documentation, include a new entry in the `NEWS.md` file.
- New functionality or bug fixes should be accompanied by a new test in the `tests/testthat` folder.
- Funcition arguments should (a) start a new line, and (b) be aligned at the equal sign, for instance:

```r
foo <- function{
  x,
  y     = 1,
  other = 2
}
```

- Updating roxygen2 comments should trigger a rebuild of the documentation using `devtools::document()`.
