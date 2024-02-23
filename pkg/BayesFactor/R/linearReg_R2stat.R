##' Using the classical R^2 test statistic for (linear) regression designs, this
##' function computes the corresponding Bayes factor test.
##'
##' This function can be used to compute the Bayes factor corresponding to a
##' multiple regression, using the classical R^2 (coefficient of determination)
##' statistic.  It can be used when you don't have access to the full data set
##' for analysis by \code{\link{lmBF}}, but you do have the test statistic.
##'
##' For details about the model, see the help for \code{\link{regressionBF}},
##' and the references therein.
##'
##' The Bayes factor is computed via Gaussian quadrature.
##' @title Use R^2 statistic to compute Bayes factor for regression designs
##' @param N number of observations
##' @param p number of predictors in model, excluding intercept
##' @param R2 proportion of variance accounted for by the predictors, excluding
##'   intercept
##' @param rscale numeric prior scale
##' @param simple if \code{TRUE}, return only the Bayes factor
##' @return If \code{simple} is \code{TRUE}, returns the Bayes factor (against the
##' intercept-only null). If \code{FALSE}, the function returns a
##' vector of length 3 containing the computed log(e) Bayes factor,
##' along with a proportional error estimate on the Bayes factor and the method used to compute it.
##' @author Richard D. Morey (\email{richarddmorey@@gmail.com}) and Jeffrey N.
##'   Rouder (\email{rouderj@@missouri.edu})
##' @keywords htest
##' @export
##' @references Liang, F. and Paulo, R. and Molina, G. and Clyde, M. A. and
##'   Berger, J. O. (2008). Mixtures of g-priors for Bayesian Variable
##'   Selection. Journal of the American Statistical Association, 103, pp.
##'   410-423
##'
##'   Rouder, J. N.  and Morey, R. D. (in press, Multivariate Behavioral Research). Bayesian testing in
##'   regression.
##'
##' @seealso \code{\link{integrate}}, \code{\link{lm}}; see
##'   \code{\link{lmBF}} for the intended interface to this function, using
##'   the full data set.
##' @examples
##' ## Use attitude data set
##' data(attitude)
##' ## Scatterplot
##' lm1 = lm(rating~complaints,data=attitude)
##' plot(attitude$complaints,attitude$rating)
##' abline(lm1)
##' ## Traditional analysis
##' ## p value is highly significant
##' summary(lm1)
##'
##' ## Bayes factor
##' ## The Bayes factor is over 400,000;
##' ## the data strongly favor hypothesis that
##' ## the slope is not 0.
##' result = linearReg.R2stat(30,1,0.6813)
##' exp(result[['bf']])

linearReg.R2stat=function(N,p,R2,rscale="medium", simple = FALSE) {
  rscale = rpriorValues("regression",,rscale)

  if(p<1)
    stop("Number of predictors must be >0")

  if(p>=(N-1))
    stop("Number of predictors must be less than N - 1 (number of data points minus 1).")

  if( (R2>=1) | (R2<0) )
    stop("Illegal R2 value (must be 0 <= R2 < 1)")

  ### Compute approximation to posterior mode of g
  ### Liang et al Eq. A.3, assuming a=b=0
  g3 = -(1 - R2) * (p + 3) #* g^3
  g2 = (N - p - 4 - 2 * (1 - R2)) #* g^2
  g1 = (N * (2 - R2) - 3) #*g
  g0 = N

  sol = polyroot(c(g0, g1, g2, g3))
  ## Pick the real solution
  modeg = Re(sol[which.min(Im(sol)^2)])
  if(modeg<=0) modeg = N/20
  log.const = integrand.regression.u(0, N, p , R2, rscaleSqr=rscale^2, log=TRUE, shift=log(modeg))

  h=integrate(integrand.regression.u,lower=-Inf,upper=Inf,N=N,p=p,R2=R2,rscaleSqr=rscale^2,log.const=log.const,shift=log(modeg))
  properror = exp(log(h[[2]]) - log(h[[1]]))
  bf = log(h$value) + log.const
  if(simple){
    return(c(B10=exp(bf)))
  }else{
    return(list(bf=bf, properror=properror, method="quadrature"))
  }
}
