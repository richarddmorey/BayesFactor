The `BayesFactor` package enables the computation of Bayes factors in standard designs, such as one- and two- sample designs, ANOVA designs, and regression.

### Installing

To install the latest stable version from CRAN, use `install.packages`:

```R
install.packages('BayesFactor')
```
or your R graphical user interface's install packages menu.

To install the latest development version, you can use `install_github` from the `devtools` package:

```R
## install devtools if necessary
install.packages('devtools')
## get BayesFactor from github
install_github('BayesFactor','richarddmorey')
```

Under Linux, you'll need standard build tools installed (`gcc`, etc).
Under OSX, you'll need Xcode with the command-line tools.
Under Windows, you'll probably have to have `Rtools` (available from CRAN).
