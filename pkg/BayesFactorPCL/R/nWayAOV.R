

##' Computes a single Bayes factor, or samples from the posterior, for an ANOVA 
##' model defined by a design matrix
##' 
##' This function is not meant to be called by end-users, although 
##' technically-minded users can call this function for flexibility beyond what 
##' the other functions in this package provide. See \code{\link{lmBF}} for a 
##' user-friendly front-end to this function. Details about the priors can be 
##' found in the help for \code{\link{anovaBF}} and the references therein.
##' 
##' Arguments \code{struc} and \code{gMap} provide a way of grouping columns of 
##' the design matrix as a factor; the effects in each group will share a common
##' \eqn{g} parameter. Only one of these arguments is needed; if both are given,
##' \code{gMap} takes precedence.
##' 
##' \code{gMap} should be a vector of the same length as the number of 
##' nonconstant rows in \code{X}. It will contain all integers from 0 to 
##' \eqn{N_g-1}{Ng-1}, where \eqn{N_g}{Ng} is the total number of \eqn{g} 
##' parameters. Each element of \code{gMap} specifies the group to which that 
##' column belongs.
##' 
##' If all columns belonging to a group are adjacent, \code{struc} can instead 
##' be used to compactly represent the groupings. \code{struc} is a vector of 
##' length \eqn{N_g}{Ng}. Each element specifies the number columns in the 
##' group. \code{gMap} is thus the \code{\link{inverse.rle}} of \code{struc}, 
##' minus 1.
##' 
##' The vector \code{rscale} should be of length \eqn{N_g}{Ng}, and contain the 
##' prior scales of the standardized effects. See Rouder et al. (2012) for more 
##' details and the help for \code{\link{anovaBF}} for some typical values.
##' 
##' The method used to estimate the Bayes factor depends on the \code{method} 
##' argument. "simple" is most accurate for small to moderate sample sizes, and 
##' uses the Monte Carlo sampling method described in Rouder et al. (2012).
##' "importance" uses an importance sampling algorithm with an importance
##' distribution that is multivariate normal on log(g). "laplace" does not
##' sample, but uses a Laplace approximation to the integral. It is expected to
##' be more accurate for large sample sizes, where MC sampling is slow. 
##' integration, and the posterior is sampled with a Gibbs sampler.
##' @title Use ANOVA design matrix to compute Bayes factors or sample posterior
##' @param y vector of observations
##' @param X design matrix whose number of rows match \code{length(y)}.
##' @param struc vector grouping the columns of \code{X} (see Details).
##' @param gMap alternative way of grouping the columns of \code{X}
##' @param rscale a vector of prior scale(s) of appropriate length (see 
##'   Details).
##' @param iterations Number of Monte Carlo samples used to estimate Bayes 
##'   factor or posterior
##' @param progress  if \code{TRUE}, show progress with a text progress bar
##' @param gibi interface for a future graphical user interface (not intended 
##'   for use by end users)
##' @param gibbs if \code{TRUE}, return samples from the posterior instead of a 
##'   Bayes factor
##' @param method the integration method (only valid if \code{gibbs=TRUE}); one 
##'   of "simple", "importance", "laplace"
##' @return If \code{posterior} is \code{FALSE}, a vector of length 2 containing
##'   the computed log(e) Bayes factor (against the intercept-only null), along 
##'   with a proportional error estimate on the Bayes factor. Otherwise, an 
##'   object of class \code{mcmc}, containing MCMC samples from the posterior is
##'   returned.
##' @export
##' @keywords htest
##' @author Richard D. Morey (\email{richarddmorey@@gmail.com}), Jeffery N. 
##'   Rouder (\email{rouderj@@missouri.edu})
##' @seealso  See \code{\link{lmBF}} for the user-friendly front end to this 
##'   function; see \code{\link{regressionBF}} and \code{anovaBF} for testing 
##'   many regression or ANOVA models simultaneously.
##' @references Rouder, J. N., Morey, R. D., Speckman, P. L., Province, J. M., 
##'   (2012) Default Bayes Factors for ANOVA Designs. Journal of Mathematical 
##'   Psychology.  56.  p. 356-374.
##' @examples
##' ## Classical example, taken from t.test() example
##' ## Student's sleep data
##' data(sleep)
##' plot(extra ~ group, data = sleep)
##' 
##' ## traditional ANOVA gives a p value of 0.00283
##' summary(aov(extra ~ group + Error(ID/group), data = sleep))
##' 
##' ## Build design matrix
##' group.column <- rep(1/c(-sqrt(2),sqrt(2)),each=10)
##' subject.matrix <- model.matrix(~sleep$ID - 1,data=sleep$ID)
##' ## Note that we include no constant column
##' X <- cbind(group.column, subject.matrix)
##' 
##' ## (log) Bayes factor of full model against grand-mean only model
##' bf.full <- nWayAOV(y = sleep$extra, X = X, struc = c(1,10), rscale=c(.5,1))
##' exp(bf.full['bf'])
##' 
##' ## Compare with lmBF result (should be about the same, give or take 1%)
##' bf.full2 <- lmBF(extra ~ group + ID, data = sleep, whichRandom = "ID")
##' bf.full2

nWayAOV<- function(y, X, struc = NULL, gMap = NULL, rscale, iterations = 10000, progress = TRUE, gibi = NULL, gibbs = FALSE, method="simple")
{
  if(!is.numeric(y)) stop("y must be numeric.")  
  if(!is.numeric(X)) stop("X must be numeric.")  
  
  N = length(y)
	X = matrix( X, nrow=N )
  
  constantCols = apply(X,2,function(v) length(unique(v)))==1
  if( sum(constantCols) > 0 )
	{
    X = X[,-which(constantCols)]
	  warning( sum(constantCols)," constant columns removed from X." )
  }
	P = ncol(X)
  
  if(!is.null(gMap)){
    if(length(gMap) != P)
      stop("Invalid gMap argument. length(gMap) must be the the same as the number of parameters (excluding intercept): ",sum(gMap)," != ",P)
    if( !all(0:max(gMap) %in% unique(gMap)) )
      stop("Invalid gMap argument: no index can be skipped.")
    
    nGs = as.integer(max(gMap) + 1)
  }else if(!is.null(struc)){
    struc = unlist(struc)
    if(sum(struc) != P)
      stop("Invalid struc argument. sum(struc) must be the the same as the number of parameters (excluding intercept): ",sum(struc)," != ",P)
    
    nGs = length(struc)
    gMap = as.integer(inverse.rle(list(values = (1:nGs)-1, lengths = struc)))
  }else{
    stop("One of gMap or struc must be defined.")
  }
  
  C = diag(N) - matrix(1/N,N,N)
  Cy = y - mean(y)
  CX = C %*% X
  XtCX = t(X) %*% CX
  XtCy = t(X) %*% C %*% as.matrix(y,cols=1)
  ytCy = var(y)*(N-1)
  
  a = rep(0.5,nGs)
  if(length(rscale)==nGs){
    b = rscale^2/2
  }else{
    stop("Length of rscale vector wrong. Was ", length(rscale), " and should be ", nGs,".")
  }
  
	nullLike = - ((N-1)/2)*log((N-1)*var(y))
	
  
  ####### Progress bar stuff
	if(!is.null(gibi)) {
		progress=TRUE;
		if(!is.function(gibi))
			stop("Malformed GIBI argument (not a function). You should not set this argument if running oneWayAOV.Gibbs from the console.")
	}
	if(progress & is.null(gibi)){
		pb = txtProgressBar(min = 0, max = 100, style = 3) 
	}else{ 
		pb=NULL 
	}
	
  pbFun = function(samps){ 
    if(progress){
    	percent = as.integer(round(samps / iterations * 100))
    	if(is.null(gibi)){
    		setTxtProgressBar(pb, percent)
    	}else{
    		gibi(percent)
    	}
    }
  }
  ############ End progress bar stuff
  
  
  if(gibbs){ # Create chains
    Z = cbind(1,X)
    ZtZ = t(Z)%*%Z
    Zty = t(Z)%*%matrix(y,ncol=1)
    
    chains = .Call("RGibbsNwayAov", 
                   as.integer(iterations), y, Z, ZtZ, Zty, as.integer(N), 
                   as.integer(P), as.integer(nGs), as.integer(gMap), rscale,
                   as.integer(iterations/100*progress), pbFun, new.env(), 
                   package="BayesFactor")
    
    dim(chains) <- c(2 + P + nGs, iterations)
    chains = mcmc(t(chains))  
    labels = c("mu",paste("beta",1:P,sep="_"),"sig2",paste("g",1:nGs,sep="_"))
    colnames(chains) = labels
    retVal = chains
 
  }else{ # Compute Bayes factor
    if(method=="simple"){
      returnList = .Call("RjeffSamplerNwayAov", 
                     as.integer(iterations), XtCX, XtCy, ytCy, 
                     as.integer(N), 
                     as.integer(P), 
                     as.integer(nGs), 
                     as.integer(gMap), a, b,
                     as.integer(iterations/100*progress), 
                     pbFun, new.env(), 
                     package="BayesFactor")
    }else if(method=="importance"){
      apx = gaussianApproxAOV(y,X,rscale,gMap)
      returnList = .Call("RimportanceSamplerNwayAov", 
                         as.integer(iterations), XtCX, XtCy, ytCy, 
                         as.integer(N), 
                         as.integer(P), 
                         as.integer(nGs), 
                         as.integer(gMap), a, b, apx$mu, apx$sig,
                         as.integer(iterations/100*progress), 
                         pbFun, new.env(), 
                         package="BayesFactor")
    }else if(method=="laplace"){
      bf = laplaceAOV(y,X,rscale,gMap)
      properror=NA
      retVal = c(bf = bf, properror=properror)
      if(inherits(pb,"txtProgressBar")) close(pb);
      return(retVal)
    }else{  
      stop("Unknown method specified.")
    }
    bf = returnList[[1]] - nullLike
    
    # estimate error
    bfSamp = returnList[[2]] - nullLike
    properror = propErrorEst(bfSamp)
    
    retVal = c(bf = bf, properror=properror)
  }
  
  if(inherits(pb,"txtProgressBar")) close(pb);
  return(retVal)
}
