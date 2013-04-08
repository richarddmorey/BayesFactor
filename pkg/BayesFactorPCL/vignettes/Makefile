# Makefile to use knitr for package vignettes

HTMLS= manual.html

all: $(HTMLS) 

clean:
	rm -rf figure *.md *.bbl *.blg *.aux *.out *.log *.png *.tex *.Rnw

%.html: %.Rmd
	$(R_HOME)/bin/Rscript -e "library(methods);library(knitr); library(BayesFactor); knit2html('$*.Rmd')"
