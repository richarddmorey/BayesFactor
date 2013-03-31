## S4 method
#####

# setMethod("plot", "BFBayesFactor", function(x, include1 = TRUE, addDenom = FALSE, sortbf=TRUE, logbase = c("log10", "log2","ln"), marginExpand=.4,pars=NULL, ...){
#   plot.BFBayesFactor(x, include1 = include1,
#                      addDenom = addDenom, 
#                      sortbf = sortbf,
#                      logbase = logbase,
#                      marginExpand = marginExpand,
#                      pars = pars, ...)
#   invisible(NULL)
# })

## S3 method
#####


#' Plot a Bayes factor object
#' 
#' This function creates a barplot of the (log) Bayes factors in a Bayes factor 
#' object. Error bars are added (though in many cases they may be too small to
#' see) in red to show the error in estimation of the Bayes factor. If a red question mark 
#' appears next to a bar, then that Bayes factor has no error estimate available.
#' @title Plot a Bayes factor object
#' @param x a BFBayesFactor object
#' @param include1 if \code{TRUE}, ensure that Bayes factor = 1 is on the plot
#' @param addDenom if \code{TRUE}, add the denominator model into the group
#' @param sortbf sort the Bayes factors before plotting them? Defaults to 
#'   \code{TRUE}
#' @param logbase the base of the log Bayes factors in the plot
#' @param marginExpand an expansion factor for the left margin, in case more 
#'   space is needed for model names
#' @param pars a list of par() settings
#' @param ... additional arguments to pass to barplot()
#' @method plot BFBayesFactor
#' @author Richard D. Morey (\email{richarddmorey@@gmail.com})
#' @examples
#' data(puzzles)
#' 
#' bfs = anovaBF(RT ~ shape*color + ID, data = puzzles, whichRandom="ID", progress=FALSE)
#' plot(bfs)
plot.BFBayesFactor <- function(x, include1=TRUE, addDenom = FALSE, sortbf=TRUE, logbase = c("log10", "log2","ln"), marginExpand = .4, pars=NULL, ...){

  # eliminate NAs
  x = x[!is.na(x)]
  
  oldPar <- par()
  on.exit(par(oldPar[c("mfrow","las",names(pars))]))
  textLogBase = logbase[1]

  logBase <- switch(textLogBase, 
                  log10=10,
                  ln=exp(1),
                  log2=2, 
                  stop('Invalid logarithm base.'))

  # Add denominator
  if(addDenom) x =  c(x, (1/x[1]) / (1/x[1]))
  if(sortbf) x = sort(x)
  
  bfs <- extractBF(x, logbf = TRUE)
  
  # Estimate left margin
  maxChar = max(nchar(rownames(bfs)))
  leftMargin = marginExpand * maxChar + 4
  
  # Errors
  whichNA = is.na(bfs$error)
  bfs$error[whichNA] = 0
  errs <- exp(bfs$bf + log(bfs$error))
  errs <- log(outer(errs,c(-1,1),'*') + exp(bfs$bf))/log(logBase)
  
  if(include1){
    rng <- range(c(0,errs))
  }else{
    rng <- range(errs)
  }
  yaxes <- seq(floor(rng[1]), ceiling(rng[2]), 1)
  ygrids <- seq(yaxes[1], yaxes[length(yaxes)], .1)

  if(textLogBase=="ln"){
    tickLab <- paste("exp(",yaxes,")",sep="")
    tickLab[yaxes==0] = "1"
  }else{
    tickLab <- logBase^yaxes
    tickLab[yaxes<0] = paste("1/",logBase^abs(yaxes[yaxes<0]),sep="")
  }


  cols = c("wheat","lightslateblue")[(bfs$bf>0) + 1]
  pars = c(pars, list(oma=c(5,leftMargin,0,1),las=1,mar=c(0,0,2,0)))
  par(pars)
  yloc <- barplot( bfs$bf/log(logBase), 
           names.arg=rownames(bfs), 
           horiz=TRUE, 
           axes=FALSE, 
           xlim=range(yaxes), 
           main = paste("vs.",x@denominator@longName), col=cols,...)

  # add error bars
  segments(errs[,1],yloc,errs[,2],yloc,col="red")
  
  # add unknown errors
  if(any(whichNA)) 
    mapply(function(x,y,adj)
      text(x,y,"?",col="red",adj=adj)
           , x=errs[whichNA,1],y=yloc[whichNA],adj=1-(errs[whichNA,1]>0))
  
  axis(1, at = yaxes, labels=tickLab, las=2)
  if(length(ygrids) < 50) abline(v=ygrids,col="gray",lty=2)
  
  abline(v=yaxes, col="gray")
  abline(v=0)
  
  invisible(NULL)
}

