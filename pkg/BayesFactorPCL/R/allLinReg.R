linearReg2 <- function(modNum, y, covariates, extraInfo = FALSE, rscale = 1, logbf = FALSE, ...){
  if(modNum==0){
    if(extraInfo){
      ret <- as.data.frame(matrix(c(1,0,0),ncol=1))
      colnames(ret) <- "null"
      return(ret)
    }else{
      return(c(null=1))
    }
  }

  y <- y - mean(y)
  covariates <- apply(covariates,2,function(v) v - mean(v))
  
  newX <- as.matrix(covariates[,binary(modNum,dim = ncol(covariates))$dicotomy])
  my.names <- colnames(covariates)[binary(modNum, dim = ncol(covariates))$dicotomy]
  
  tX <- t(newX)
  beta <- solve(tX%*%newX)%*%tX%*%y
  resid <- y - newX%*%beta
  
  R2 <- 1 - sum(resid^2)/sum(y^2)  
  p <- ncol(newX)
  N <- nrow(newX)
  
  bf <- linearReg.Quad(N=N, p=p, R2=R2, rscale, logbf, ...)
  label <- paste(my.names,collapse=' + ')
  
  if(extraInfo){
    ret <- as.data.frame(matrix(c(bf,R2,p),ncol=1))
    colnames(ret) <- label
    return(ret)
  }else{
    names(bf) <- label
    return(bf)
  }  
}

allLinearReg <- function(y, covariates, extraInfo = FALSE, progress = FALSE, rscale = 1, logbf = FALSE, ...){
  topModel <- 2 ^ ncol(covariates) - 1
  ret <- list()
  progressEvery <- round((topModel + 1) / 100)
  if(progressEvery < 1) progressEvery <- 1
  if(identical(progress, TRUE)){
    pb <- txtProgressBar(min=0, max=100, style=3) 
  }
  for(mod in 0:topModel){
    ret <- c(ret, linearReg2(mod, y, covariates, extraInfo, rscale, logbf, ...))
    if(!identical(progress, FALSE)){
      if( ( (mod + 1)%%progressEvery ) == 0){
        if(identical(progress, TRUE)){
          setTxtProgressBar(pb,round(mod/topModel*100))
        }else if(is.function(progress)){
          progress( round(100 * mod / topModel) )
        }
      }
    }
  }
  if(identical(progress, TRUE)){
    if(inherits(pb,"txtProgressBar")) close(pb)
  }  
  if(extraInfo){
    my.names <- names(ret)
    ret <- matrix(unlist(ret), nrow = 3)
    ret <- data.frame(t(ret))
    rownames(ret) <- my.names
    colnames(ret) <- c("bf","R2","p")
    return(ret[order(ret$bf, decreasing=TRUE),])
  }else{
    return(sort(unlist(ret),decreasing=TRUE))
  }
}



