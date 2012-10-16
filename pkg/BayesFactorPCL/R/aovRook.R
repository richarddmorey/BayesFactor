rookEnv <- new.env(parent=emptyenv())
  
aovGUI <- function(y,dataFixed=NULL,dataRandom=NULL){
  # clean up old GUI
  if(exists("aov$s",envir=rookEnv)){
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
  rookEnv$aov$bfEnv$doneBFs = c(null=0)
  rookEnv$aov$bfEnv$nFac = dim(dataFixed)[2]
  rookEnv$aov$bfEnv$designMatrices = list()
  rookEnv$aov$bfEnv$designMatrices[[ 2 ^ rookEnv$aov$bfEnv$nFac ]] = matrix(nrow=0,ncol=0)
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

setJSONdata <- function(req, res){
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
  if(req$GET()$what=="analysis"){
    if(is.null(req$GET()$model)){
      res$write(toJSON("Error"))  
      res$finish()
      return()
    }else{
      modNum <- as.integer(req$GET()$model)
    }   
    
       if( modNum==0 & is.null(rookEnv$aov$bfEnv$dataRandom)){
           returnList <- list(
             model = 0,
             name = "null",
             bf = 1,
             isBase = FALSE,
             iterations = "",
             rscaleFixed = "",
             rscaleRandom = "",
             duration = "",
             time = format(Sys.time(), "%H:%M:%S"),
             timestamp = as.integer(Sys.time())
           )
         res$write(toJSON(
           returnList
         ))  
       }
       
       rscaleFixed <- ifelse (is.null(req$GET()$rscaleFixed), 0.5, as.numeric(req$GET()$rscaleFixed))
       rscaleRandom <- ifelse (is.null(req$GET()$rscaleRandom), 1, as.numeric(req$GET()$rscaleRandom))
       iterations <- ifelse (is.null(req$GET()$iterations), 10000, as.integer(req$GET()$iterations))
       duration <- system.time({
          bf <- nWayAOV2(modNum, env = rookEnv$aov$bfEnv, 
               rscaleFixed = rscaleFixed, rscaleRandom = rscaleRandom, 
               iterations = iterations)[1]
       })[[3]]
       modelName = ifelse(modNum==0, "null", names(bf))
       
       rookEnv$aov$bfEnv$doneBFs[[modelName]] = as.numeric(bf)
       
       returnList <- list(
            model = modNum,
            name = modelName,
            bf = as.numeric(bf),
            isBase = FALSE,
            iterations = iterations,
            rscaleFixed = rscaleFixed,
            rscaleRandom = rscaleRandom,
            duration = duration,
            time = format(Sys.time(), "%H:%M:%S"),
            timestamp = as.integer(Sys.time())
         )
        res$write(toJSON(
            returnList
       ))  
  }
  res$finish()
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
      }else{
        setJSONdata(req, res)
      }
    },
    '^/.*\\.png$' = function(env){
      req <- Request$new(env)
      res <- Response$new()
      res$header('Content-type','image/png')
      if (length(rookEnv$aov$bfEnv$doneBFs) == 0){
        res$finish()
        return()
      } 
      baseBF <- ifelse(is.null(req$GET()$baseBF), 0, as.numeric(req$GET()$baseBF))
      bfs <- unlist(rookEnv$aov$bfEnv$doneBFs) - baseBF
      
      t <- tempfile()
      png(file=t)
      png(t,width=800,height=300)
      
      rng <- range(bfs/log(10))
      yaxes <- seq(floor(rng[1]), ceiling(rng[2]), 1)
      ygrids <- seq(yaxes[1], yaxes[length(yaxes)], .1)
      
      par(mar=c(4,20,1,1),las=1)
      barplot(sort(bfs/log(10)),horiz=TRUE, axes=FALSE, xlim=range(yaxes))
      axis(1, at = yaxes, lab=10^yaxes)
      abline(v=0)
      abline(v=ygrids,col="gray",lty=2)
      dev.off()
      
      
      res$body <- t
      names(res$body) <- 'file'
      res$finish()
    },
  '.*' = Redirect$new('/www/aov.html')
  )
)