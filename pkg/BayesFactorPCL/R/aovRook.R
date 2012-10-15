rookEnv <- new.env(parent=emptyenv())
  
aovGUI <- function(y,dataFixed=NULL,dataRandom=NULL){
  # clean up old GUI
  if(exists("aov",envir=rookEnv)){
    stopAovGUI()
  }
  
  rookEnv$aov <- new.env(parent = rookEnv)
  
  # Convert to factors if needed.
  dataFixed = data.frame(dataFixed)
  if(any(!sapply(dataFixed,is.factor)) & !is.null(dataFixed)){
    dataFixed <- data.frame(lapply(dataFixed, as.factor))
    warning("Converted columns of dataFixed to factors.")
  } 
  dataRandom = data.frame(dataRandom)
  if(any(!sapply(dataRandom,is.factor)) & !is.null(dataRandom)){
    dataRandom <- data.frame(lapply(dataRandom, as.factor))
    warning("Converted columns of dataRandom to factors.")
  } 
  
  rookEnv$aov$bfEnv = new.env(parent = rookEnv$aov)
  rookEnv$aov$bfEnv$nFac = dim(dataFixed)[2]
  rookEnv$aov$bfEnv$designMatrices = list()
  rookEnv$aov$bfEnv$dataFixed = dataFixed
  rookEnv$aov$bfEnv$y = y
  rookEnv$aov$bfEnv$totalN = length(as.vector(y))
  rookEnv$aov$bfEnv$dataRandom = dataRandom  
  
  rookEnv$aov$bfEnv$allEffects = sapply(1:(2^rookEnv$aov$bfEnv$nFac-1),
                                        other.design,env=rookEnv$aov$bfEnv,type='n')
  
  rookEnv$aov$s <- Rhttpd$new()
  ## Not run: 
  rookEnv$aov$s$start(quiet=TRUE)
  rookEnv$aov$s$add(name="aov",
        app=aovApp)
  rookEnv$aov$s$browse("aov")
}

stopAovGUI <- function(){
  rookEnv$aov$s$remove(all=TRUE)
  rm(aov, envir=rookEnv)
}

aovApp <- Builder$new(
  Static$new(
    urls = '/www',
    root = system.file('.', package='BayesFactor')
  ),
  URLMap$new(
    '^/data' = function(env){
      req <- Request$new(env)
      res <- Response$new()
      if (is.null(req$GET()$what)){
        res$finish()
        return()
      }
      if(req$GET()$what=="fixed"){
        res$write(toJSON(
          rookEnv$aov$bfEnv$allEffects
            ))  
      }
      if(req$GET()$what=="random"){
        res$write(toJSON(
          names(rookEnv$aov$bfEnv$dataRandom)
        ))  
      }
      if(req$GET()$what=="nFac"){
        res$write(toJSON(
          rookEnv$aov$bfEnv$nFac
        ))  
      }
      res$finish()
    },
    '^/.*\\.png$' = function(env){
      req <- Request$new(env)
      res <- Response$new()
      res$header('Content-type','image/png')
      if (is.null(req$GET()$n)){
        n <- 100
      } else {
        n <- as.integer(req$GET()$n)
      }
      t <- tempfile()
      png(file=t)
      png(t,width=200,height=200)
      par(mar=rep(0,4))
      plot(rnorm(n),col=rainbow(n,alpha=runif(n,0,1)),pch='.',cex=c(2,3,4,5,10,50))
      dev.off()
      res$body <- t
      names(res$body) <- 'file'
      res$finish()
    },
  '.*' = Redirect$new('/www/aov.html')
  )
)