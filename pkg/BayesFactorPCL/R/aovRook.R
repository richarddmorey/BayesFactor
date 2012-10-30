rookEnv <- new.env(parent=emptyenv())
RservePort <- 6401
RserveArgs <- "--no-save"

getURL <- function(loc){
  con = url(loc, open = "r")
  readLines(con)
  close(con)
}

RserveCleanup <- function(){
  if(rookEnv$RserveStatus=="free" & exists("RserveSession", envir=rookEnv)){
    rookEnv$RserveConnection = RSattach(rookEnv$RserveSession)
    rm("RserveSession", envir=rookEnv)
    RSshutdown(rookEnv$RserveConnection)
    RSclose(rookEnv$RserveConnection)
    rm("RserveConnection", envir=rookEnv)
  }
}

RserveAov2 <- function(tokens,updateURL,...){
  tryStart = try( { Rserve(port=RservePort,args=RserveArgs) }, silent=TRUE )
  if(inherits(tryStart, "try-error")) {
    stopAovGUI()
    stop("Could not start Rserve.")
    return()
  }
  rookEnv$RserveConnection <- RSconnect(port=RservePort)
  
  theseArgs = list(...)
  RSassign(rookEnv$RserveConnection,theseArgs)
  RSassign(rookEnv$RserveConnection,tokens)
  RSassign(rookEnv$RserveConnection,updateURL)
  RSeval(rookEnv$RserveConnection, quote(library(BayesFactor, quietly=TRUE)))

  cmd = parse(text='
    lastUpdate = 0
    for(i in 1:length(tokens)){
      token <- names(tokens)[i]
      model <- tokens[i]
      gibi <- function(percent){
        rightNow = proc.time()[3]
        timeSince = rightNow - lastUpdate
        if( timeSince > 0.5){
          myURL = paste(updateURL, "?token=",token,"&percent=",percent,"&status=running&since=",timeSince, sep="")
          BayesFactor:::getURL(myURL)
          lastUpdate <<- rightNow
        }
      }
      theseArgs$modNum = model
      theseArgs$progress = TRUE
      theseArgs$gibi=gibi
      duration <- system.time({
        bf <- do.call("nWayAOV2", theseArgs)[1]
      })[[3]]
      myURL = paste(updateURL, "?token=",token,
                    "&percent=100&status=done&bf=",bf,
                    "&duration=",as.numeric(duration),
                    "&time=",as.integer(Sys.time()),sep="")
      BayesFactor:::getURL(myURL)
  }
  myURL = paste(updateURL, "?status=alldone",sep="")
  BayesFactor:::getURL(myURL)
  ')
  RSassign(rookEnv$RserveConnection,cmd)
  rookEnv$RserveStatus = "running"
  rookEnv$RserveSession <- RSevalDetach(rookEnv$RserveConnection, "eval(cmd)")
}

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

aovGUI <- function(y,dataFixed=NULL,dataRandom=NULL, iterations = 1000, rscaleFixed=.5, rscaleRandom=1){
  if(!exists("aov",envir=rookEnv)){
    rookEnv$aov <- new.env(parent = rookEnv)
  }
  
  rookEnv$RserveStatus = "free"
  
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
  
  rookEnv$aov$defaults = list(iterations=iterations, rscaleRandom=rscaleRandom, rscaleFixed=rscaleFixed)
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
  if(exists("RserveConnection",envir=rookEnv)){
    RSshutdown(rookEnv$RserveConnection)
    RSclose(rookEnv$RserveConnection)
    rm("RserveConnection",envir=rookEnv)
  }
  
  if(exists("s",envir=rookEnv$aov)){
    rookEnv$aov$s$remove(all=TRUE)
    rm("aov", envir=rookEnv)
  }
}

setJSONdata <- function(req, res){
  if(req$GET()$what=="fixed"){
    res$write(toJSON(
      list(effects=rookEnv$aov$bfEnv$allEffects,
           nFac=rookEnv$aov$bfEnv$nFac
      )
    ))
    return()
  }
  if(req$GET()$what=="random"){
    res$write(toJSON(
      names(rookEnv$aov$bfEnv$dataRandom)
    ))  
    return()
  }
  if(req$GET()$what=="nFac"){
    res$write(toJSON(
      rookEnv$aov$bfEnv$nFac
    ))
    return()
  }
  if(req$GET()$what=="analysis"){
    if(rookEnv$RserveStatus!="free"){
      res$write(toJSON(
        list(status="busy")
      )) 
      return()
    }
    
    tokens <- unlist(strsplit(req$params()$tokens, ","))
    models <- unlist(strsplit(req$params()$models, ","))
    
    if(is.null(tokens)){
      res$write(toJSON(
        list(token=-1)
      )) 
      return()
    }  
    tokensExist = sapply(tokens,function(e){!is.null(rookEnv$aov$status[[e]])})
    tokenContent <- rookEnv$aov$status[tokens]
    
    if(all(tokensExist)){
      res$write(toJSON(
        list(status="oldtokens")
      ))
      return()
    } 
    
    if(is.null(models)){
      res$write(toJSON(
        list(token=-1)
      ))  
      return()
    }  
    
    tokens = tokens[!tokensExist]
    models = models[!tokensExist]
    
    rscaleFixed <- ifelse (is.null(req$GET()$rscaleFixed), rookEnv$aov$defaults$rscaleFixed, as.numeric(req$GET()$rscaleFixed))
    rscaleRandom <- ifelse (is.null(req$GET()$rscaleRandom), rookEnv$aov$defaults$rscaleRandom, as.numeric(req$GET()$rscaleRandom))
    iterations <- ifelse (is.null(req$GET()$iterations), rookEnv$aov$defaults$iterations, as.integer(req$GET()$iterations))
    
    tokensToAnalyze = c()
    modelsToAnalyze = c()
    
    for(i in 1:length(models)){
      modNum <- as.integer(models[i])
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
          token = tokens[i],
          status = "done"
        )
        rookEnv$aov$status[[ tokens[i] ]] <- newToken(tokens[i], returnList)
        rookEnv$aov$status[[ tokens[i] ]]$status = "done"
        rookEnv$aov$status[[ tokens[i] ]]$percent = 100
        rookEnv$aov$status[[ tokens[i] ]]$finish = as.integer(Sys.time())
      }else{
        tokensToAnalyze = c(tokensToAnalyze, tokens[i])
        modelsToAnalyze = c(modelsToAnalyze, modNum)
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
          token = tokens[i]
        )
        rookEnv$aov$status[[ tokens[i] ]] <- newToken(tokens[i], returnList)
      }
    }        
      updateURL = paste(rookEnv$aov$s$full_url(1),"/rserve",sep="")
      names(modelsToAnalyze) = tokensToAnalyze
      #print(modelsToAnalyze)
    
      RserveAov2(modelsToAnalyze,updateURL, env=rookEnv$aov$bfEnv,
                 iterations = iterations,
                 rscaleFixed = rscaleFixed,
                 rscaleRandom = rscaleRandom)
         
    res$write(toJSON(list(status="started")))
    return()
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
      RserveCleanup()
      if (is.null(req$GET()$what)){
        return(res$finish())
      }else{
        setJSONdata(req, res)
        return(res$finish())
      }
    },
    '^/update' = function(env){
      req <- Request$new(env)
      res <- Response$new()
      RserveCleanup()
      #cat("Update request:",req$query_string(),"\n")
      if (is.null(req$GET()$tokens)){
        res$write(toJSON(-1))
        return(res$finish())
      }else{  
        tokens <- unlist(strsplit(req$params()$tokens, ","))
        percents <- sapply(tokens, function(e){ rookEnv$aov$status[[e]]$percent},USE.NAMES = FALSE)
        statuses <- sapply(tokens, function(e){ rookEnv$aov$status[[e]]$status},USE.NAMES = FALSE)
        allstatus <- ifelse(all(statuses=="done"),"done","running")
        returnLists <- sapply(tokens, function(e){ toJSON(rookEnv$aov$status[[e]]$returnList)},USE.NAMES = FALSE)      
        retJSON = list(
          tokens=tokens, percents=percents, statuses=statuses, allstatus=allstatus, returnLists=returnLists
        )
        res$write(toJSON(retJSON,asIs=TRUE))
        return(res$finish())
      }
    },
    '^/rserve' = function(env){
      req <- Request$new(env)
      res <- Response$new()
      RserveCleanup()
      #cat("Rserve request:",req$query_string(),"\n")
      status = req$GET()$status
      percent = req$GET()$percent
      token = req$GET()$token
      if(is.null(token) & status=="alldone"){
        rookEnv$RserveStatus = "free"
      }
      if(!is.null(token) & status=="done"){
        bf = req$GET()$bf
        finish = req$GET()$time
        duration = req$GET()$duration
        rookEnv$aov$status[[token]]$returnList$bf <- as.numeric(bf)
        rookEnv$aov$status[[token]]$returnList$duration <- duration
        rookEnv$aov$status[[token]]$status = "done"
        rookEnv$aov$status[[token]]$percent = 100
        rookEnv$aov$status[[token]]$finish = as.integer(finish)        
      }
      if(!is.null(token) & status=="running"){
        rookEnv$aov$status[[token]]$status = "running"
        rookEnv$aov$status[[token]]$percent = req$GET()$percent
      }
      res$write(0)
      return(res$finish())
    },
    '^/getdefaults' = function(env){
      req <- Request$new(env)
      res <- Response$new()
      RserveCleanup()
      res$write(toJSON(rookEnv$aov$defaults))
      res$finish()
    },
    '^/saveobj' = function(env){
      req <- Request$new(env)
      res <- Response$new()
      RserveCleanup()
            
      if (is.null(req$params()$BFobj)){
        return(res$finish())
      }
      textLogBase <- ifelse(is.null(req$params()$logBase), "log10", req$params()$logBase)
      logBase <- switch(textLogBase, log10=10,ln=exp(1),log2=2)
      

      # Parse Bayes factors from JSON and put them in data.frame
      bfs <- fromJSON(req$GET()$BFobj)
      if(length(bfs)>1){
        bfs <- merge_recurse(lapply(bfs,data.frame))
      }else{
        bfs <- data.frame(bfs) 
      }
      
      baseBF <- bfs$bf[bfs$isBase]
      bfs[[paste("logbf.",textLogBase,sep="")]] <- (bfs$bf - baseBF) / log(logBase)
      bfs$bf <- exp(bfs$bf - baseBF)
      
      if(req$params()$which=="CSV"){
        t <- tempfile()
        write.csv(bfs,file=t)
        res$header('Content-type','text/csv')
        res$header("Content-Disposition", "attachment;filename=bfs.csv")
        res$body <- t
        names(res$body) <- 'file'
        return(res$finish())
      }else if(req$params()$which=="R"){
        .GlobalEnv$BayesFactorTable = bfs
        return(res$finish())
      }else{
        return(res$finish())
      }
    },
    '^/bfs.png' = function(env){
      req <- Request$new(env)
      res <- Response$new()
      RserveCleanup()
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
      bfs = bfs[!is.na(bfs$bf),]
      
      
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
