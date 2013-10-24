cd ..
rm ABCoptim_*
R CMD REMOVE ABCoptim
R CMD build ABCoptim
R CMD INSTALL ABCoptim_0.13*
R CMD check --as-cran ABCoptim_*
cd ABCoptim
