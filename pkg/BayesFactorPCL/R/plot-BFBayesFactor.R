## S4 method
#####

setMethod("plot", "BFBayesFactor", function(x, include1 = TRUE, addDenom = FALSE, sortbf=TRUE, logbase = c("log10", "log2","ln"), marginExpand=.4,pars=NULL, ...){
  plot.BFBayesFactor(x, include1 = include1,
                     addDenom = addDenom, 
                     sortbf = sortbf,
                     logbase = logbase,
                     marginExpand = marginExpand,
                     pars = pars, ...)
  invisible(NULL)
})

## S3 method
#####

plot.BFBayesFactor <- function(x, include1=TRUE, addDenom = FALSE, sortbf=TRUE, logbase = c("log10", "log2","ln"), marginExpand = .4, pars=NULL, ...){

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
  pars = c(pars, list(mar=c(4,leftMargin,4,1),las=1))
  par(pars)
  yloc <- barplot( bfs$bf/log(logBase), 
           names.arg=rownames(bfs), 
           horiz=TRUE, 
           axes=FALSE, 
           xlim=range(yaxes), 
           main = paste("vs.",x@denominator@longName), col=cols,...)

  # add error bars
  segments(errs[,1],yloc,errs[,2],yloc,col="red")
  
  axis(1, at = yaxes, labels=tickLab, las=2)
  if(length(ygrids) < 50) abline(v=ygrids,col="gray",lty=2)
  
  abline(v=yaxes, col="gray")
  abline(v=0)
  
  invisible(NULL)
}

