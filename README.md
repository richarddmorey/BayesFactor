`BayesFactor` is an R package that enables the computation of Bayes factors in standard designs, such as one- and two- sample designs, ANOVA designs, and regression.

### Installing

To install the latest stable version from CRAN, use `install.packages`:

```R
install.packages('BayesFactor', dependencies = TRUE)
```
or your R graphical user interface's install packages menu.

To install the latest development version, you can use `install_github` from the `devtools` package:

```R
## install devtools if necessary
install.packages('devtools')
## Load devtools package for install_github()
library(devtools)
## get BayesFactor from github
install_github('richarddmorey/BayesFactor', subdir='pkg/BayesFactor', dependencies = TRUE)
```

Under Linux, you'll need standard build tools installed (`gcc`, etc).
Under OSX, you'll need Xcode with the command-line tools.
Under Windows, you'll probably have to have `Rtools` (available from CRAN).
