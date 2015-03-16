

#'Functions to compute Bayes factor hypothesis tests for common research designs
#'and hypotheses.
#'
#'This package contains function to compute Bayes factors for a number of 
#'research designs and hypotheses, including t tests, ANOVA, and linear 
#'regression, and contingency tables.
#'
#'\tabular{ll}{ Package: \tab BayesFactor\cr Type: \tab Package\cr Version: \tab
#'0.9.11\cr Date: \tab 2015-3-16\cr License: \tab GPL 2.0\cr LazyLoad: \tab 
#'yes\cr } The following methods are currently implemented, with more to follow:
#'
#'general linear models (including linear mixed effects models): \code{\link{generalTestBF}}, \code{\link{lmBF}}
#'
#'linear regression: \code{\link{regressionBF}}, \code{\link{lmBF}}, 
#'\code{\link{linearReg.R2stat}};
#'
#'t tests: \code{\link{ttestBF}}, \code{\link{ttest.tstat}};
#'
#'meta-analytic t tests: \code{\link{meta.ttestBF}}
#'
#'ANOVA: \code{\link{anovaBF}}, \code{\link{lmBF}}, \code{\link{oneWayAOV.Fstat}};
#'
#'contingency tables: \code{\link{contingencyTableBF}};
#'
#'single proportions: \code{\link{proportionBF}};
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
#'@author Richard D. Morey and Jeffrey N. Rouder (with contributions from Tahira Jamil)
#'  
#'  Maintainer: Richard D. Morey <richarddmorey@@gmail.com>
#'@references Liang, F. and Paulo, R. and Molina, G. and Clyde, M. A. and 
#'  Berger, J. O. (2008). Mixtures of g-priors for Bayesian Variable Selection. 
#'  Journal of the American Statistical Association, 103, pp. 410-423
#'  
#'  Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., \& Iverson, G. 
#'  (2009). Bayesian t-tests for accepting and rejecting the null hypothesis. 
#'  Psychonomic Bulletin & Review, 16, 225-237
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

#'Hraba and Grant (1970) children's doll preference data
#'
#'Hraba and Grant (1970) describe a replication of Clark and Clark (1947) in which
#'black and white children from Lincoln, Nebraska were shown dolls that were either black
#'or white. They were then asked a series of questions, including "Give me the doll that is 
#'a nice doll." This data set contains the frequency of children giving the same-race or different race doll in
#'response to this question.
#'@name raceDolls
#'@docType data
#'@format A matrix with 2 rows and 2 columns. Rows give doll preference; colums give the 
#'race of the child.
#'@source Hraba, J. and Grant, G. (1970). Black is Beautiful: A reexamination of 
#'racial preference and identification. Journal of Personality and Social Psychology, 16, 398-402.
#'
#'@keywords datasets
#'@examples
#'
#'data(raceDolls)
#'
#'## chi-square test
#'## Barely significant with continuity correction
#'chisq.test(raceDolls)
#'
#'## Bayes factor test (assuming independent binomial sampling plan)
#'## Very little evidence for the alternative of lack of independence
#'bf = contingencyTableBF(raceDolls, sampleType = "indepMulti", fixedMargin = "cols")
#'bf
NULL
