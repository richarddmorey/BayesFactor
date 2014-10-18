
##' Using the classical F test statistic for a balanced one-way design, this function computes the corresponding Bayes factor test.
##' 
##' For F statistics computed from balanced one-way designs, this function can
##' be used to compute the Bayes factor testing the model that all group means
##' are not equal to the grand mean, versus the null model that all group means
##' are equal. It can be used when you don't have access to the full data set
##' for analysis by \code{\link{lmBF}}, but you do have the test statistic.
##' 
##' For details about the model, see the help for \code{\link{anovaBF}}, and the references therein.
##' 
##' The Bayes factor is computed via Gaussian quadrature.
##' @title Use F statistic to compute Bayes factor for balanced one-way designs
##' @param F F statistic from classical ANOVA
##' @param N number of observations per cell or group
##' @param J number of cells or groups
##' @param rscale numeric prior scale
##' @param simple if \code{TRUE}, return only the Bayes factor
##' @return If \code{simple} is \code{TRUE}, returns the Bayes factor (against the 
##' intercept-only null). If \code{FALSE}, the function returns a 
##' vector of length 2 containing the computed log(e) Bayes factor,
##' along with a proportional error estimate on the Bayes factor.
##' @export
##' @keywords htest
##' @author Richard D. Morey (\email{richarddmorey@@gmail.com})
##' @references Morey, R. D., Rouder, J. N., Pratte, M. S., \& Speckman, P. L.
##'   (2011). Using MCMC chain outputs to efficiently estimate Bayes factors.
##'   Journal of Mathematical Psychology, 55, 368-378
##'   
##' @note \code{oneWayAOV.Fstat} should only be used with F values obtained from
##'   balanced designs.
##' @examples
##' ## Example data "InsectSprays" - see ?InsectSprays
##' require(stats); require(graphics)
##' boxplot(count ~ spray, data = InsectSprays, xlab = "Type of spray", 
##'         ylab = "Insect count", main = "InsectSprays data", varwidth = TRUE, 
##'         col = "lightgray")
##' 
##' ## Classical analysis (with transformation)
##' classical <- aov(sqrt(count) ~ spray, data = InsectSprays)
##' plot(classical)
##' summary(classical)
##' 
##' ## Bayes factor (a very large number)
##' Fvalue <- anova(classical)$"F value"[1]
##' result <- oneWayAOV.Fstat(Fvalue, N=12, J=6)
##' exp(result[['bf']])
##' @seealso \code{\link{integrate}}, \code{\link{aov}}; see \code{\link{lmBF}} for the intended interface to this function, using the full data set.

oneWayAOV.Fstat = function(F, N, J, rscale="medium", simple = FALSE)
{
  rscale = rpriorValues("allNways","fixed",rscale)
  log.const = marginal.g.oneWay(1,F=F,N=N,J=J,rscale=rscale,log=TRUE)
  integral = integrate(marginal.g.oneWay,lower=0,upper=Inf,F=F,N=N,J=J,rscale=rscale,log.const=log.const)
  properror = exp(log(integral[[2]]) - log(integral[[1]]))
	bf = log(integral[[1]]) + log.const
	if(simple){
	  return(c(B10=exp(bf)))
  }else{
    return(c(bf=bf, properror=properror))
	}
}
