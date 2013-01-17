\name{plot.BFBayesFactorTop}
\alias{plot.BFBayesFactorTop}
\title{Plot a Bayes factor top-down object}
\usage{
  \method{plot}{BFBayesFactorTop} (x, include1 = TRUE,
    addDenom = FALSE, sortbf = FALSE,
    logbase = c("log10", "log2", "ln"), marginExpand = 0.4,
    pars = NULL, ...)
}
\arguments{
  \item{x}{a BFBayesFactorTop object}

  \item{include1}{if \code{TRUE}, ensure that Bayes factor
  = 1 is on the plot}

  \item{addDenom}{if \code{TRUE}, add the denominator model
  into the group}

  \item{sortbf}{sort the Bayes factors before plotting
  them? Defaults to \code{TRUE}}

  \item{logbase}{the base of the log Bayes factors in the
  plot}

  \item{marginExpand}{an expansion factor for the left
  margin, in case more space is needed for model names}

  \item{pars}{a list of par() settings}

  \item{...}{additional arguments to pass to barplot()}
}
\description{
  Plot a Bayes factor top-down object
}
\details{
  This function creates a barplot of the (log) Bayes
  factors in a Bayes factor object. Error bars are added
  (though in many cases they may be too small to see) in
  red to show the error in estimation of the Bayes factor.
}
\examples{
data(puzzles)

bfs = anovaBF(RT ~ shape*color + ID, data = puzzles, whichRandom="ID", whichModels='top', progress=FALSE)
plot(bfs)
}
\author{
  Richard D. Morey (\email{richarddmorey@gmail.com})
}

