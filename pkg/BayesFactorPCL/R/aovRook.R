rookEnv <- new.env(parent=emptyenv())
  
newToken <- function(token, returnList)
{
  newToken <- list(status="running",
                   percent=0,
                   start=as.integer(Sys.time()),
                   finish=NULL,
                   returnList=returnList,
                   token=token)
  return(newToken)
}

aovGUI <- function(y,dataFixed=NULL,dataRandom=NULL){
  if(!exists("aov",envir=rookEnv)){
    rookEnv$aov <- new.env(parent = rookEnv)
  }
  
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
  rookEnv$aov$bfEnv$designMatrices[[ 2 ^ rookEnv$aov$bfEnv$nFac ]] = matrix(nrow=0,ncol=0)
  rookEnv$aov$bfEnv$dataFixed = dataFixed
  rookEnv$aov$bfEnv$y = y
  rookEnv$aov$bfEnv$totalN = length(as.vector(y))
  rookEnv$aov$bfEnv$dataRandom = dataRandom  
  rookEnv$aov$status <- list()
  
  rookEnv$aov$bfEnv$allEffects = sapply(1:(2^rookEnv$aov$bfEnv$nFac-1),
                                        other.design,env=rookEnv$aov$bfEnv,type='n')
  
  if(!exists("s",envir=rookEnv$aov)){
    rookEnv$aov$s <- Rhttpd$new()
    rookEnv$aov$s$start(quiet=TRUE)
    rookEnv$aov$s$add(name="aov",
          app=aovApp)
  }
  rookEnv$aov$s$browse("aov")
}

stopAovGUI <- function(){
  rookEnv$aov$s$remove(all=TRUE)
  rm(aov, envir=rookEnv)
}

setJSONdata <- function(req, res){
  if(req$GET()$what=="fixed"){
    res$write(toJSON(
      list(effects=rookEnv$aov$bfEnv$allEffects,
           nFac=rookEnv$aov$bfEnv$nFac
      )
    ))
    return(res$finish())
  }
  if(req$GET()$what=="random"){
    res$write(toJSON(
      names(rookEnv$aov$bfEnv$dataRandom)
    ))  
    return(res$finish())
  }
  if(req$GET()$what=="nFac"){
    res$write(toJSON(
      rookEnv$aov$bfEnv$nFac
    ))
    return(res$finish())
  }
  if(req$GET()$what=="analysis"){
    token <- req$GET()$token
    model <- req$GET()$model
    
    if(!is.null(token)){
      tokenContent <- rookEnv$aov$status[[token]]
      if(is.null(tokenContent)){
        res$write(toJSON(
          list(token=-1)
        ))  
        return(res$finish())
      }else{
        res$write(toJSON(
          rookEnv$aov$status[[token]]
        ))
        return(res$finish())
      }
    }
    if(is.null(model)){
      res$write(toJSON(
        list(token=-1)
      ))  
      return(res$finish())
    }  
    
    token <- substring(tempfile(tmpdir=""),5)
    
    rscaleFixed <- ifelse (is.null(req$GET()$rscaleFixed), 0.5, as.numeric(req$GET()$rscaleFixed))
    rscaleRandom <- ifelse (is.null(req$GET()$rscaleRandom), 1, as.numeric(req$GET()$rscaleRandom))
    iterations <- ifelse (is.null(req$GET()$iterations), 10000, as.integer(req$GET()$iterations))
    
    
    modNum <- as.integer(model)
    modelName = ifelse(modNum==0, "null", 
                       paste(joined.design(modNum, env=rookEnv$aov$bfEnv, other="n"), collapse=" + "))
    
    if( modNum==0 & is.null(rookEnv$aov$bfEnv$dataRandom)){
      returnList <- list(
        model = 0,
        name = "null",
        niceName = "null",
        bf = 1,
        isBase = FALSE,
        iterations = "",
        rscaleFixed = "",
        rscaleRandom = "",
        duration = "",
        time = format(Sys.time(), "%H:%M:%S"),
        timestamp = as.integer(Sys.time()),
        token = token,
        status = "done"
      )
      rookEnv$aov$status[[token]]$status = "done"
      rookEnv$aov$status[[token]]$percent = 100
      rookEnv$aov$status[[token]]$finish = as.integer(Sys.time())
    }else{
      returnList <- list(
        model = modNum,
        name = modelName,
        niceName = modelName,
        bf = NULL,
        isBase = FALSE,
        iterations = iterations,
        rscaleFixed = rscaleFixed,
        rscaleRandom = rscaleRandom,
        duration = "",
        time = format(Sys.time(), "%H:%M:%S"),
        timestamp = as.integer(Sys.time()),
        token = token
      )
    }
    
    gibi <- function(percent){
      rookEnv$aov$status[[token]]$percent <- percent
    }
    
    rookEnv$aov$status[[token]] <- newToken(token, returnList)
    res$write(toJSON(
      rookEnv$aov$status[[token]]
    ))
    
    finish <- res$finish()
    if(rookEnv$aov$status[[token]]$status == "done") return(finish)
    
    
    duration <- system.time({
      bf <- nWayAOV2(modNum, env = rookEnv$aov$bfEnv, 
                     rscaleFixed = rscaleFixed, rscaleRandom = rscaleRandom, 
                     iterations = iterations, progress=TRUE, gibi=gibi)[1]
    })[[3]]
    
    
    rookEnv$aov$status[[token]]$returnList$bf <- as.numeric(bf)
    rookEnv$aov$status[[token]]$returnList$duration <- duration
    rookEnv$aov$status[[token]]$status = "done"
    rookEnv$aov$status[[token]]$percent = 100
    rookEnv$aov$status[[token]]$finish = as.integer(Sys.time())
    
    return(finish)
  }
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
        return(res$finish())
      }else{
        finish <- setJSONdata(req, res)
        return(finish)
      }
    },
    '^/.*\\.png$' = function(env){
      req <- Request$new(env)
      res <- Response$new()
      res$header('Content-type','image/png')
      if (is.null(req$params()$BFobj)){
        return(res$finish())
      }
      textLogBase <- ifelse(is.null(req$params()$logBase), "log10", req$params()$logBase)
      logBase <- switch(textLogBase, log10=10,ln=exp(1),log2=2)
      
      # Parse Bayes factors from JSON and put them in data.frame
      bfs <- fromJSON(req$GET()$BFobj)
      bfs <- merge_recurse(lapply(bfs,data.frame))
      
      baseBF <- bfs$bf[bfs$isBase]
      bfs$bf <- bfs$bf - baseBF
      
      t <- tempfile()
      png(file=t)
      png(t,width=700,height=350)
      
      rng <- range(bfs$bf/log(logBase))
      yaxes <- seq(floor(rng[1]), ceiling(rng[2]), 1)
      ygrids <- seq(yaxes[1], yaxes[length(yaxes)], .1)
      
      if(textLogBase=="ln"){
        tickLab <- paste("exp(",yaxes,")",sep="")
        tickLab[yaxes==0] = "1"
      }else{
        tickLab <- logBase^yaxes
        tickLab[yaxes<0] = paste("1/",logBase^abs(yaxes[yaxes<0]),sep="")
      }
      
      par(mar=c(4,20,1,1),las=1)
      barplot(sort(bfs$bf/log(logBase)), names.arg=bfs$niceName[order(bfs$bf)], horiz=TRUE, axes=FALSE, xlim=range(yaxes))
      axis(1, at = yaxes, lab=tickLab, las=2)
      abline(v=0)
      if(length(ygrids) < 50) abline(v=ygrids,col="gray",lty=2)
      dev.off()
      
      
      res$body <- t
      names(res$body) <- 'file'
      res$finish()
    },
  '.*' = Redirect$new('/www/aov.html')
  )
)