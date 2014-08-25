oneWayAOV.Gibbs = function(y,iterations=10000,rscale="medium", progress=TRUE, callback=function(...) as.integer(0), logbf=FALSE){
  
  rscale = rpriorValues("allNways","fixed",rscale)
  N = as.integer(colSums(!is.na(y)))
  J=as.integer(dim(y)[2])
  I=as.integer(dim(y)[1])
  iterations = as.integer(iterations)
  if(progress){
    pb = txtProgressBar(min = 0, max = 100, style = 3) 
  }else{ 
    pb=NULL 
  }
  
  
  pbFun = function(samps){ 
    if(progress){
      percent = as.integer(round(samps / iterations * 100))
      setTxtProgressBar(pb, percent)
    }
  }
  
  output = .Call("RgibbsOneWayAnova", y, N, J, I, rscale, iterations,
                 progress, pbFun, callback, new.env(), package="BayesFactor")
  
  if(progress) close(pb);
  rownames(output[[1]]) = c("mu",paste("beta",1:J,sep=""),"CMDESingle","CMDEDouble","sig2","g")			
  names(output[[2]])=c("logCMDESingle","logCMDEDouble","logCMDESingleKahan","logCMDEDoubleKahan")
  
  logPriorDensDouble = dmvnorm(rep(0,J),rep(0,J),diag(J),log=TRUE)  
  
  logPostDensDouble = logMeanExpLogs(output[[1]][1+J+2,])
  lbf = logPostDensDouble - logPriorDensDouble
  
  if(logbf){
    return(list(chains=mcmc(t(output[[1]])), BF=-lbf))
  }else{
    return(list(chains=mcmc(t(output[[1]])), BF=exp(-lbf)))
  }
  
  
}



marginal.g.oneWay = function(g,F,N,J,rscale)
{
  dfs = (J-1)/(N*J-J)
  omega = (1+(N*g/(dfs*F+1)))/(N*g+1)
  m = log(rscale) - 0.5*log(2*pi) - 1.5*log(g) - rscale^2/(2*g) - (J-1)/2*log(N*g+1) - (N*J-1)/2*log(omega)
  exp(m)
}

