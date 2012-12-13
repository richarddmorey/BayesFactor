\name{allNways}
\alias{allNways}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the Bayes factors in an ANOVA design
}
\description{
This function computes Bayes factors for all main-effects and interaction contrasts in an ANOVA design. 
}
\usage{
allNways(y, dataFixed = NULL, dataRandom = NULL, iterations = 10000, 
        whichModels = "withmain", progress = TRUE, rscaleFixed = .5, 
        rscaleRandom = 1, logbf = FALSE, multicore = FALSE, 
        extraInfo=FALSE, only.top=NULL, ...) 
    
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	A vector of observations	
}
  \item{dataFixed}{
A data frame containing the fixed factors, one column per factor. Columns will be converted to factors if they are not factors already.
}
  \item{dataRandom}{
A data frame containing the random factors, one column per factor. Columns will be converted to factors if they are not factors already.
}
  \item{iterations}{
Number of Monte Carlo samples to perform (see details).  
}
  \item{whichModels}{
One of three values: 'all', 'top', or 'withmain' (see details).
}
\item{progress}{
If \code{TRUE}, a text progress bar is displayed. Currently unavailable if \code{multicore=TRUE}.
}
\item{rscaleFixed}{The prior scale of effect size for factors that are treated as a fixed.
}
\item{rscaleRandom}{
The prior scale of effect size for factors that are treated as random.}

\item{logbf}{
If \code{TRUE}, return the natural logarithm of the Bayes factor
}
\item{multicore}{
If \code{TRUE}, and the \code{\link[doMC:doMC-package]{doMC}} 
 package is installed, the function will attempt to use multiple cores. Number of cores is set by options(cores).
}
\item{extraInfo}{
  If \code{TRUE}, return a data frame containing useful information. If \code{FALSE}, only the Bayes factors are returned. 
}
\item{only.top}{
This argument does nothing; it has been replaced by \code{whichModels}, and will be removed in the future.
}

\item{...}{
Additional arguments to pass to \code{\link{nWayAOV.MC}}
}
}

\details{
The models, priors, and methods of computation are provided in Rouder et al. (2012).  	
The ANOVA model for a vector of observations \eqn{y} is
\deqn{ y = \mu + X_1 \theta_1 + \ldots + X_p\theta_p +\epsilon,}
where \eqn{\theta_1,\ldots,\theta_p} are vectors of main-effect and interaction effects, \eqn{X_1,\ldots,X_p} are corresponding design matrices, and \eqn{\epsilon} is a vector of zero-centered noise terms with variance \eqn{\sigma^2}.  Zellner and Siow (1980) inspired g-priors are placed on effects, but with a separate g-prior parameter for each covariate:
\deqn{\theta_1~N(0,g_1\sigma^2), \ldots,  \theta_p~N(0,g_p \sigma^2).}  A Jeffries prior is placed on \eqn{\mu} and \eqn{\sigma^2}.  Independent scaled inverse-chi-square priors with one degree of freedom are placed on \eqn{g_1,\ldots,g_p}.  The square-root of the scale for g's corresponding to fixed and random effects is given by \code{rscaleFixed} and \code{rscaleRandom}, respectively.

When a factor is treated as random, there are as many main effect terms in the vector \eqn{\theta} as levels.  When a factor is treated as fixed, the sums-to-zero linear constraint is enforced by centering the corresponding design matrix, and there is one fewer main effect terms as levels.  The Cornfield-Tukey model of interactions is assumed.  Details are provided in Rouder et al. (2012)

Bayes factors are computed by integrating the likelihood with respect to the priors on parameters.  The integration of all parameters except \eqn{g_1,\ldots,g_p} may be expressed in closed-form; the integration of \eqn{g_1,\ldots,g_p} is performed through Monte Carlo sampling, and \code{iterations} is the number of iterations used to estimate the Bayes factor.  

\code{allNways} computes Bayes factors for either all submodels or select submodels missing a single main effect or covariate, depending on the argument \code{only.top}.  Bayes factors for specific models are computed by \code{\link{nWayAOV.MC}}, an options to this back-end function may be passed.  See \code{\link{nWayAOV.MC}} for more details.

If no random factors are specified, the null model assumed by \code{allNways} is the grand-mean only model. If random factors are specified, the null model is the model with all random factors plus a grand mean. Thus, \code{allNways} cannot currently be used to test random factors. Testing random factors is possible with \code{\link{nWayAOV.MC}}.

The argument \code{whichModels} controls which models are tested. Setting \code{whichModels} to 'all' will test all models that can be created by including or not including a main effect or interaction. 'top' will test all models that can be created by removing or leaving in a main effect or interaction term from the full model. 'withmain' will test all models, with the constraint that if an interaction is included, the corresponding main effects are also included.
}


\value{
If \code{whichModels} is 'all' or 'withmain', a vector of Bayes factors, each measuring the support for a model relative to the null model. 

If \code{whichModels} is 'top', a vector of Bayes factors, each measuring the support for a model relative to the full model. Only the models that can be formed by removing an effect or interaction from the full model are analyzed and returned.
}

\references{
Rouder, J. N., Morey, R. D., Speckman, P. L., Province, J. M., (2012) Default Bayes Factors for ANOVA Designs. Journal of Mathematical Psychology.  56.  p. 356-374.

Zellner, A. and Siow, A., (1980) Posterior Odds Ratios for Selected Regression Hypotheses.  In Bayesian Statistics: Proceedings of the First Interanational Meeting held in Valencia (Spain).  Bernardo, J. M., Lindley, D. V., and Smith A. F. M. (eds), pp. 585-603.  University of Valencia.

}
\author{
Richard D. Morey, Jeffrey N. Rouder
}



%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{
The function \code{allNways} will compute Bayes factors for all possible combinations of fixed factors and interactions, against the null hypothesis that \emph{all} effects are 0. The total number of tests computed will be \eqn{2^{2^K - 1}}{2^(2^K - 1)} for \eqn{K} fixed factors. This number increases very quickly with the number of factors. For instance, for a five-way ANOVA, the total number of tests exceeds two billion. Even though each test takes a fraction of a second, the time taken for all tests could exceed your lifetime.

It is possible to reduce the number of models tested by only testing the most complex model and every restriction that can be formed by removing one factor or interaction using the \code{whichModels} argument. Setting this argument to 'top' reduces the number of tests to \eqn{2^K-1}, which is more manageable. The Bayes factor for each restriction against the most complex model can be interpreted as a test of the removed factor/interaction. Setting \code{whichModels} to 'withmain' will not reduce the number of tests as much as 'top' but the results may be more interpretable, since an interaction is only allowed when all interacting effects (main or interaction) are also included in the model.
}
\seealso{
\code{\link{nWayAOV.MC}}, \code{\link{aov}}, \code{\link{ttest.Quad}}, \code{\link{oneWayAOV.Quad}},\code{\link{oneWayAOV.Gibbs}}, \code{\link[doMC:doMC-package]{doMC}}, \code{\link{options}}
}
\examples{
## Classical example, taken from t.test() example
## Student's sleep data
data(sleep)
plot(extra ~ group, data = sleep)

## traditional ANOVA gives a p value of 0.00283
summary(aov(extra ~ group + Error(ID/group), data = sleep))

## Gives a Bayes factor of about 11.6
## in favor of the alternative hypothesis
allNways(y = sleep$extra,
      dataFixed = data.frame(sleep$group), 	   
      dataRandom = data.frame(sleep$ID), progress=FALSE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
