##' Computes a single Bayes factor, or samples from the posterior, for an ANOVA 
##' model defined by a design matrix
##' 
##' This function is not meant to be called by end-users, although 
##' technically-minded users can call this function for flexibility beyond what 
##' the other functions in this package provide. See \code{\link{lmBF}} for a 
##' user-friendly front-end to this function. Details about the priors can be 
##' found in the help for \code{\link{anovaBF}} and the references therein.
##' 
##' Argument \code{gMap} provides a way of grouping columns of 
##' the design matrix as a factor; the effects in each group will share a common
##' \eqn{g} parameter. \code{gMap} should be a vector of the same length as the number of 
##' nonconstant rows in \code{X}. It will contain all integers from 0 to 
##' \eqn{N_g-1}{Ng-1}, where \eqn{N_g}{Ng} is the total number of \eqn{g} 
##' parameters. Each element of \code{gMap} specifies the group to which that 
##' column belongs.
##' 
##' If all columns belonging to a group are adjacent, \code{struc} can instead 
##' be used to compactly represent the groupings. \code{struc} is a vector of 
##' length \eqn{N_g}{Ng}. Each element specifies the number columns in the 
##' group. 
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
##' be more accurate for large sample sizes, where MC sampling is slow. If 
##' \code{method="auto"}, then an initial run with both samplers is done, and 
##' the sampling method that yields the least-variable samples is chosen. The 
##' number of initial test iterations is determined by 
##' \code{options(BFpretestIterations)}.
##' 
##' If posterior samples are requested, the posterior is sampled with a Gibbs 
##' sampler.
##' @title Use ANOVA design matrix to compute Bayes factors or sample posterior
##' @param y vector of observations
##' @param X design matrix whose number of rows match \code{length(y)}.
##' @param gMap vector grouping the columns of \code{X} (see Details).
##' @param rscale a vector of prior scale(s) of appropriate length (see 
##'   Details).
##' @param iterations Number of Monte Carlo samples used to estimate Bayes 
##'   factor or posterior
##' @param progress  if \code{TRUE}, show progress with a text progress bar
##' @param callback callback function for third-party interfaces 
##' @param gibbs will be deprecated. See \code{posterior}
##' @param posterior if \code{TRUE}, return samples from the posterior using 
##' Gibbs sampling, instead of the Bayes factor
##' @param ignoreCols if \code{NULL} and \code{posterior=TRUE}, all parameter
##'   estimates are returned in the MCMC object. If not \code{NULL}, a vector of
##'   length P-1 (where P is number of columns in the design matrix) giving which
##'   effect estimates to ignore in output
##' @param thin MCMC chain to every \code{thin} iterations. Default of 1 means 
##'   no thinning. Only used if \code{posterior=TRUE}
##' @param method the integration method (only valid if \code{posterior=FALSE}); one 
##'   of "simple", "importance", "laplace", or "auto"
##' @param continuous either FALSE if no continuous covariates are included, or 
##'   a logical vector of length equal to number of columns of X indicating 
##'   which columns of the design matrix represent continuous covariates
##' @param noSample if \code{TRUE}, do not sample, instead returning NA. This is 
##'   intended to be used with functions generating and testing many models at one time, 
##'   such as \code{\link{anovaBF}}
##' @return If \code{posterior} is \code{FALSE}, a vector of length 2 containing
##'   the computed log(e) Bayes factor (against the intercept-only null), along 
##'   with a proportional error estimate on the Bayes factor. Otherwise, an 
##'   object of class \code{mcmc}, containing MCMC samples from the posterior is
##'   returned.
##' @note Argument \code{struc} has been deprecated. Use \code{gMap}, which is the \code{\link{inverse.rle}} of \code{struc}, 
##' minus 1.
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
##' bf.full <- nWayAOV(y = sleep$extra, X = X, gMap = c(0,rep(1,10)), rscale=c(.5,1))
##' exp(bf.full[['bf']])
##' 
##' ## Compare with lmBF result (should be about the same, give or take 1%)
##' bf.full2 <- lmBF(extra ~ group + ID, data = sleep, whichRandom = "ID")
##' bf.full2

nWayAOV<- function(y, X, gMap, rscale, iterations = 10000, progress = options()$BFprogress, callback = function(...) as.integer(0), gibbs = NULL, posterior = FALSE, ignoreCols=NULL, thin=1, method="auto", continuous=FALSE, noSample = FALSE)
{  
  if(!is.numeric(y)) stop("y must be numeric.")  
  if(!is.numeric(X)) stop("X must be numeric.")  
  if(!is.function(callback)) stop("Invalid callback.")  
  
  if(!is.null(gibbs)){
    warning("Argument 'gibbs' to nWayAOV will soon be deprecated. Use 'posterior' instead.")
    posterior = gibbs
  }
     
  N = length(y)
	X = matrix( X, nrow=N )
  
  constantCols = apply(X,2,function(v) length(unique(v)))==1 & !apply(X,2,function(v) all(v==0))
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
  }else{
    stop("gMap must be defined.")
  }
  
  if(is.null(ignoreCols)) ignoreCols = rep(0,P)
  
  if(length(rscale)!=nGs){
    stop("Length of rscale vector wrong. Was ", length(rscale), " and should be ", nGs,".")
  }
      
  # Rearrange design matrix if continuous columns are included
  # We will undo this later if chains have to be returned
  if(!identical(continuous,FALSE)){
    if(all(continuous)){
      #### If all covariates are continuous, we want to use Gaussian quadrature.
        Cy = matrix(y - mean(y), ncol=1)
        CX = t(t(X) - colMeans(X))
        R2 = t(Cy)%*%CX%*%solve(t(CX)%*%CX)%*%t(CX)%*%Cy / (t(Cy)%*%Cy)
        bf = linearReg.R2stat(N=N,p=ncol(CX),R2=R2,rscale=rscale)  
        return(bf)
    }
    if(length(continuous) != P) stop("argument continuous must have same length as number of predictors")
    if(length(unique(gMap[continuous]))!=1) stop("gMap for continuous predictors don't all point to same g value")
    
    # Sort chains so that continuous covariates are together, and first
    sortX = order(!continuous)
    revSortX = order(sortX)
    X = X[,sortX]
    gMap = gMap[sortX]
    ignoreCols = ignoreCols[sortX]
    continuous = continuous[sortX]
    incCont = sum(continuous)
  }else{
    revSortX = sortX = 1:ncol(X)
    incCont = as.integer(0)
  }
  
  # What if we can use quadrature?
  if(nGs==1 & !posterior & all(!continuous)) 
    return(singleGBayesFactor(y,X,rscale,gMap, incCont))
  
  if(posterior)
    return(nWayAOV.Gibbs(y, X, gMap, rscale, iterations, incCont, sortX, revSortX, progress, ignoreCols, thin, continuous, noSample, callback))
 
  if(!noSample){
    if(method %in% c("simple","importance","auto")){
      return(doNwaySampling(method, y, X, rscale, 
                   iterations, gMap, incCont, progress, callback))
    }else if(method=="laplace"){
      bf = laplaceAOV(y,X,rscale,gMap,incCont)
      return(c(bf = bf, properror=NA))
    }else{  
      stop("Unknown method specified.")
    }
  }
  
  return(c(bf = NA, properror=NA))
}


nWayAOV.Gibbs = function(y, X, gMap, rscale, iterations, incCont, sortX, revSortX, progress, ignoreCols, thin, continuous, noSample, callback)
{
  P = ncol(X)
  nGs = as.integer( max(gMap) + 1 )

  # Check thinning to make sure number is reasonable
  if( (thin<1) | (thin>(iterations/3)) ) stop("MCMC thin parameter cannot be less than 1 or greater than iterations/3. Was:", thin)
  
  if(!identical(continuous,FALSE)){
    if(all(continuous)){
        chains = linearReg.Gibbs(y, X, iterations = iterations, 
                                 rscale = rscale, progress = progress, callback = callback)
        return(chains)
    } 
    if(length(continuous) != P) stop("argument continuous must have same length as number of predictors")
    if(length(unique(gMap[continuous]))!=1) stop("gMap for continuous predictors don't all point to same g value")
  }
  
  nOutputPars = sum(1-ignoreCols)
  
  if(noSample){ # Return structure of chains
    chains = matrix(NA,2,nOutputPars + 2 + nGs)
  }else{  
    chains = jzs_Gibbs(iterations, y, cbind(1,X), rscale, 1, gMap, table(gMap), incCont, FALSE, 
                       as.integer(ignoreCols), as.integer(thin), as.logical(progress), callback, 1)
  }
  chains = mcmc(chains)  
  
  # Unsort the chains if we had continuous covariates
  if(incCont){
    # Account for ignored columns when resorting
    revSort = 1+order(sortX[!ignoreCols])
    chains[,1 + 1:nOutputPars] = chains[,revSort]
    labels = c("mu",paste("beta",1:P,sep="_")[!ignoreCols[revSortX]],"sig2",paste("g",1:nGs,sep="_"))
  }else{
    labels = c("mu",paste("beta",1:P,sep="_")[!ignoreCols],"sig2",paste("g",1:nGs,sep="_"))
  }
  colnames(chains) = labels
  return(chains)
}

