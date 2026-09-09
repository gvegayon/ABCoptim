.PHONY: build install check docs help

help:
	@echo "Available targets:"
	@echo "  build     Build the package"
	@echo "  install   Install the package"
	@echo "  check     Check the package"
	@echo "  docs      Generate package documentation"

build:
	Rscript -e 'devtools::build()'

install:
	Rscript -e 'devtools::install()'

check:
	Rscript -e 'devtools::check()'

docs:
	Rscript -e 'devtools::document()'

README.md: README.qmd
	quarto render README.qmd