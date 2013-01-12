# Makefile to use knitr for package vignettes

HTMLS= index.html

all: $(HTMLS) 

clean:
	rm -rf figure *.md *.bbl *.blg *.aux *.out *.log *.png *.tex

%.html: %.Rmd
	$(R_HOME)/bin/Rscript -e "library(knitr); library(BayesFactor); knit2html('$*.Rmd')"
