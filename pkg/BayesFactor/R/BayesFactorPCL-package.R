

#'Functions to compute Bayes factor hypothesis tests for common research designs
#'and hypotheses.
#'
#'This package contains function to compute Bayes factors for a number of 
#'research designs and hypotheses, including t tests, ANOVA, and linear 
#'regression.
#'
#'\tabular{ll}{ Package: \tab BayesFactor\cr Type: \tab Package\cr Version: \tab
#'0.9.5\cr Date: \tab 2013-3-20\cr License: \tab GPL 2.0\cr LazyLoad: \tab 
#'yes\cr } The following methods are currently implemented, with more to follow:
#'
#'Linear regression: \code{\link{regressionBF}} \code{\link{lmBF}}, 
#'\code{\link{linearReg.R2stat}};
#'
#'t test: \code{\link{ttestBF}}, \code{\link{ttest.tstat}};
#'
#'ANOVA: \code{\link{anovaBF}}, \code{\link{lmBF}}, 
#'\code{\link{oneWayAOV.Fstat}};
#'
#'Other useful functions: \code{\link{posterior}}, for sampling from posterior 
#'distributions; \code{\link{recompute}}, for re-estimating a Bayes factor or 
#'posterior distribution; \code{\link{compare}}, to compare two model
#'posteriors; and \code{\link{plot.BFBayesFactor}}, for plotting Bayes factor
#'objects.
#'
#'@name BayesFactor-package
#'@aliases BayesFactor-package BayesFactor
#'@docType package
#'@author Richard D. Morey and Jeffrey N. Rouder
#'  
#'  Maintainer: Richard D. Morey <richarddmorey@@gmail.com>
#'@seealso \code{\link[BAS:BAS-package]{BAS}}
#'@references Liang, F. and Paulo, R. and Molina, G. and Clyde, M. A. and 
#'  Berger, J. O. (2008). Mixtures of g-priors for Bayesian Variable Selection. 
#'  Journal of the American Statistical Association, 103, pp. 410-423
#'  
#'  Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., \& Iverson, G. 
#'  (2009). Bayesian t-tests for accepting and rejecting the null hypothesis. 
#'  Psychonomic Bulletin & Review, 16, 752-760
#'  
#'  Rouder, J. N., Morey, R. D., Speckman, P. L., Province, J. M., (2012) 
#'  Default Bayes Factors for ANOVA Designs. Journal of Mathematical Psychology.
#'  56.  p. 356-374.
#'  
#'  Perception and Cognition Lab (University of Missouri): Bayes factor 
#'  calculators. \url{http://pcl.missouri.edu/bayesfactor}
#'@keywords htest
#'@examples
#'
#'## See specific functions for examples.
#'
#'@useDynLib BayesFactor
NULL





#'Puzzle completion times from Hays (1994)
#'
#'Puzzle completion time example data from Hays (1994).
#'
#'Hays (1994; section 13.21, table 13.21.2, p. 570) describes a experiment
#'wherein 12 participants complete four puzzles each. Puzzles could be either
#'square or round, and either monochromatic or in color. Each participant
#'completed every combination of the two factors.
#'
#'@name puzzles
#'@docType data
#'@format A data frame with 48 observations on 3 variables.  \describe{
#'\item{RT}{Puzzle completion time, in minutes} \item{ID}{the
#'subject identifier} \item{shape}{shape of the puzzle (round or
#'square)} \item{color}{color content of the puzzle (monochromatic or
#'color)} }
#'@source Hays, W. L. (1994), Statistics (5th edition), Harcourt Brace, Fort
#'Worth, Texas
#'@keywords datasets
#'@examples
#'
#'data(puzzles)
#'
#'## classical ANOVA
#'## Both color and shape are significant, interaction is not
#'classical <- aov(RT ~ shape*color + Error(ID/(shape*color)), data=puzzles)
#'summary(classical)
#'
#'## Bayes Factor
#'## Best model is main effects model, no interaction
#' anovaBF(RT ~ shape*color + ID, data = puzzles, whichRandom = "ID", progress=FALSE)
#'
#'
NULL



