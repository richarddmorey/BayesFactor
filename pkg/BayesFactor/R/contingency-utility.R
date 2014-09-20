# By Tahira Jamil, edited by Richard Morey (July 2014)
# Bayes Factor from Gunel & Dickey Paper, For different sampling models 
# (BF1: Poisson; BF2: Multinomial; BF3: Poisson; BF4: Hypergeometric)
#The Bayes factor are provided in favour of alternative hypothesis(against the null hypothesis: no association)

##############################################
## Utility functions
##############################################
ldirich <- function(a) {
  val <- sum(lgamma(a)) - lgamma(sum(a))
  return(val)
}        


ldirich1 <- function(y, a) {
  val <- sum(lgamma(a) - lgamma(a+y))
  return(val)
}


pcomb <- function( x, log=TRUE ) {
  x <- as.vector(x)
  n <- sum(x)
  ans = lgamma(n+1) - sum(lgamma(x+1))
  if(log){
    return(ans)
  }else{
    return(exp(ans))
  }
}


#########################################################################################
# Bayes Factor Gunel & Dickey  Equation 4.2
# Poisson Sampling
#########################################################################################
contingencyPoisson<-function (y, a){
  n <- sum(y)
  d <- dim(y)
  I <- d[1]
  J <- d[2]
  b <- I*J*a/n
  a <- a + 0 * y
  ac <- colSums(a)
  ar <- rowSums(a)
  yc <- colSums(y)
  yr <- rowSums(y)
  oc <- 1 + 0 * yc
  or <- 1 + 0 * yr
  
  lbf<-lgamma(sum(a) - (I-1)*(J-1)) - ((I-1)*(J-1)*log(1 + 1/b)) -
    lgamma(sum(y) + sum(a) - (I-1)*(J-1)) - ldirich(yr + ar-(J- 1)*or) +
    ldirich( ar - (J - 1) * or) - ldirich(yc + ac - (I - 1) * oc) +
    ldirich( ac - (I - 1) * oc) - ldirich1(y,a)
  return(lbf)
}



#########################################################################################
# Bayes Factor Gunel & Dickey  Equation 4.4
# Joint Multinomial Sampling
#########################################################################################

contingencyJointMultinomial <-function (y, a){
  a <- a + 0 * y
  ac <- colSums(a)
  ar <- rowSums(a)
  yc <- colSums(y)
  yr <- rowSums(y)
  d <- dim(y)
  oc <- 1 + 0 * yc
  or <- 1 + 0 * yr
  I <- d[1]
  J <- d[2]
  lbf <- ldirich(c(y) + c(a)) + ldirich(ar - (J - 1) * or) + 
    ldirich(ac - (I - 1) * oc) - ldirich(c(a)) - 
    ldirich(yr + ar - (J - 1) * or) - ldirich(yc + ac - (I - 1) * oc)
    
  return(lbf)
}


#########################################################################################
# Bayes Factor Gunel & Dickey  Equation 4.7
#Binomial/ Independent Multinomial Sampling
#########################################################################################

contingencyIndepMultinomial<-function (y, a){
  
  a <- a + 0 * y
  ac <- colSums(a)
  ar <- rowSums(a)
  yc <- colSums(y)
  yr <- rowSums(y)
  d <- dim(y)
  oc <- 1 + 0 * yc
  or <- 1 + 0 * yr
  I <- d[1]
  J <- d[2]
  lbf <- ldirich(ac - (I - 1) * oc) + ldirich(ar) + ldirich(c(y) + c(a)) -
    ldirich(yc + ac - (I - 1) * oc) - ldirich(yr + ar) - ldirich(c(a))
  
  return(lbf)
}


#########################################################################################
# Bayes Factor Gunel & Dickey  Equation 4.11
# Hypergeometric Condition on both Margins
#########################################################################################
contingencyHypergeometric<-function (y, a) {
  
  if(!identical(dim(y),as.integer(c(2,2)))) stop("hypergeometric contingency tables restricted to 2 x 2 tables; see help for contingencyTableBF()")
  
  a <- a + 0 * y
  ac <- colSums(a)
  ar <- rowSums(a)
  yc <- colSums(y)
  yr <- rowSums(y)
  d <- dim(y)
  oc <- 1 + 0 * yc
  or <- 1 + 0 * yr
  I <- d[1]
  J <- d[2]
  
  sumg<-function(y,a) {
    M <- c(yc,yr) 
    stopifnot(M[1]+M[2] == M[3]+M[4]) # To check both marginal totals are equal
    
    upper <- min(M[c(1,3)])
    lower <- 0
    if (min(M) < upper) lower <- upper - min(M)
    all.M <- t(sapply(lower:upper, function(i) c(a=i, b=M[1] - i, c=M[3] - i,
                                                 d=M[4] - M[1] + i)))
    a <- a + 0 * y
    n.sim<-upper+1
    g1<- rep(NA,n.sim)
    
    for (s in 1:n.sim)  {
      y1<-matrix(all.M[s,],2,2)
      g1[s]<-pcomb(y1) + ldirich(c(y1)+c(a)) - ldirich(c(a))
    }
    
    return (logMeanExpLogs(g1) + log(length(g1)))
  }
  
  sum.mar<- sumg(y,a)
  lbf<-ldirich(c(y) + c(a)) + pcomb(yc) + pcomb(yr) - ldirich(c(a)) - sumg(y,a)
  
  return(lbf)  
}

###########################
## Sampling
## All code below by Richard Morey
###########################

sampleContingency <- function(model, type, fixedMargin, prior, data, iterations, ...)                  
{
  if(type == "poisson"){
    if(model == "non-independence"){
      chains = samplePoissonContingencyAlt(prior, data, iterations, ...)
    }else if(model == "independence"){
      chains = samplePoissonContingencyNull(prior, data, iterations, ...)
    }
  }else if(type == "joint multinomial"){
    if(model == "non-independence"){
      chains = sampleJointMultiContingencyAlt(prior, data, iterations, ...)
    }else if(model == "independence"){
      chains = sampleJointMultiContingencyNull(prior, data, iterations, ...)
    }
  }else if(type == "independent multinomial"){
    if(model == "non-independence"){
      chains = sampleIndepMultiContingencyAlt(fixedMargin, prior, data, iterations, ...)
    }else if(model == "independence"){
      chains = sampleIndepMultiContingencyNull(fixedMargin, prior, data, iterations, ...)
    }
  }else if(type == "hypergeometric"){
    if(model == "non-independence"){
      chains = sampleHypergeomContingencyAlt(prior, data, iterations, ...)
    }else if(model == "independence"){
      chains = sampleHypergeomContingencyNull(prior, data, iterations, ...)
    }
  }else{
    stop("Unknown model type.")
  }
}

samplePoissonContingencyNull <- function(prior, data, iterations, noSample=FALSE, ...)
{  
  a = prior
  b = length(data) * a / sum(data)
  I = nrow(data)
  J = ncol(data)

  if(noSample){
    samples = data.frame(matrix(as.numeric(NA), 1, I*J + 1 + I + J))
  }else{
    lambda = rgamma(iterations, sum(data) + I*J*(a - 1) + I + J - 1, b + 1)
    pi_i = rdirichlet(iterations, rowSums(data) + a - J + 1)
    pi_j = rdirichlet(iterations, colSums(data) + a - I + 1)
    lambda_ij = t(sapply( 1:iterations, function(i) as.vector( lambda[i] * outer( pi_i[i,], pi_j[i,] ) ) ) )
    samples = data.frame( lambda_ij, lambda, pi_i, pi_j )
  }
  
  cn1 = paste0("lambda[",outer(1:nrow(data), 1:ncol(data),paste,sep=","),"]")
  cn2 = paste0("pi[",1:nrow(data),",*]")
  cn3 = paste0("pi[*,",1:ncol(data),"]")
  colnames(samples) = c(cn1,"lambda..",cn2,cn3)
  
  return(mcmc(samples))
}

samplePoissonContingencyAlt <- function(prior, data, iterations, noSample=FALSE, ...)
{
  
  data = as.matrix(data)
  a = prior
  IJ = length(data)
  b = IJ * a / sum(data)

  a.post = data + a
  b.post = data*0 + b + 1
  
  if(noSample){
    samples = data.frame(matrix(as.numeric(NA), 1, IJ))
  }else{
    samples = rgamma(iterations * IJ, a.post, rate = b.post)
    dim(samples) = c(IJ, iterations)
    samples = data.frame(t(samples))
  }
  
  cn = paste0("lambda[",outer(1:nrow(data), 1:ncol(data),paste,sep=","),"]")
  colnames(samples) = cn
  
  return(mcmc(samples))

}


sampleJointMultiContingencyNull <- function(prior, data, iterations, noSample = FALSE, ...)
{
  a = prior
  I = nrow(data)
  J = ncol(data)
  
  if(noSample){
    samples = data.frame(matrix(as.numeric(NA), 1, I*J + I + J))
  }else{
    pi_i = rdirichlet(iterations, rowSums(data) + a - J + 1)
    pi_j = rdirichlet(iterations, colSums(data) + a - I + 1)
    pi_ij = t(sapply( 1:iterations, function(i) as.vector( outer( pi_i[i,], pi_j[i,] ) ) ) )
    samples = data.frame( pi_ij, pi_i, pi_j )
  }
  
  cn1 = paste0("pi[",outer(1:nrow(data), 1:ncol(data),paste,sep=","),"]")
  cn2 = paste0("pi[",1:nrow(data),",*]")
  cn3 = paste0("pi[*,",1:ncol(data),"]")
  colnames(samples) = c(cn1,cn2,cn3)
  
  return(mcmc(samples))
}


sampleJointMultiContingencyAlt <- function(prior, data, iterations, noSample = FALSE, ...)
{
  a = prior
  I = nrow(data)
  J = ncol(data)
  
  if(noSample){
    samples = data.frame(matrix(as.numeric(NA), 1, I*J))
  }else{
    pi_ij = rdirichlet( iterations, as.matrix(data) + a )
    samples = data.frame( pi_ij )
  }
  
  cn = paste0("pi[",outer(1:nrow(data), 1:ncol(data),paste,sep=","),"]")
  colnames(samples) = cn
  
  return(mcmc(samples))
  
}


sampleIndepMultiContingencyNull <- function(fixedMargin, prior, data, iterations, noSample = FALSE, ...)
{
  a = prior
  I = nrow(data)
  J = ncol(data)
  
  if(noSample){
    if(fixedMargin == "rows"){
      samples = data.frame(matrix(as.numeric(NA), 1, I*J + I + J))
    }else{
      samples = data.frame(matrix(as.numeric(NA), 1, I*J + J + J))
    }
  }else{
    if(fixedMargin == "rows"){
      pi_star = rdirichlet( iterations, rowSums(data) + J*a )
      omega = rdirichlet( iterations, colSums(data) + a )
      pi_ij = t(sapply(1:iterations, 
                       function(i){
                         as.vector(outer( pi_star[i,], omega[i,] ))
                       }))
    }else{
      pi_star = rdirichlet( iterations, colSums(data) + I*a )
      omega = rdirichlet( iterations, rowSums(data) + a )
      pi_ij = t(sapply(1:iterations, 
                       function(i){
                         as.vector(outer( omega[i,], pi_star[i,] ))
                       }))      
    }
    samples = data.frame( pi_ij, pi_star, omega )
  }
  
  cn1 = paste0("pi[",outer(1:nrow(data), 1:ncol(data),paste,sep=","),"]")
  if(fixedMargin == "rows"){
    cn2 = paste0("pi[",1:nrow(data),",*]")
    cn3 = paste0("omega[*,",1:ncol(data),"]")
  }else{
    cn2 = paste0("pi[*,",1:ncol(data),"]") 
    cn3 = paste0("omega[",1:nrow(data),",*]")
  }
  colnames(samples) = c(cn1,cn2,cn3)
  
  return(mcmc(samples))
}


sampleIndepMultiContingencyAlt <- function(fixedMargin, prior, data, iterations, noSample = FALSE, ...)
{
  a = prior
  I = nrow(data)
  J = ncol(data)
  
  if(noSample){
    if(fixedMargin == "rows"){
      samples = data.frame(matrix(as.numeric(NA), 1, I*J + I + I*J))
    }else{
      samples = data.frame(matrix(as.numeric(NA), 1, I*J + J + I*J))
    }
  }else{
    if(fixedMargin == "rows"){
      pi_star = rdirichlet( iterations, rowSums(data) + J*a )
      omega = t(replicate(iterations, as.vector(t(apply(data, 1, function( v ) rdirichlet( 1, v + a ))))))
      pi_ij = t(sapply(1:iterations, 
                       function(i){
                         as.vector(rep(pi_star[i,], J) * omega[i,])
                       }))
    }else{
      pi_star = rdirichlet( iterations, colSums(data) + I*a )
      omega = t(replicate(iterations, as.vector(apply(data, 2, function( v ) rdirichlet( 1, v + a )))))
      pi_ij = t(sapply(1:iterations, 
                       function(i){
                         as.vector(rep(pi_star[i,], each = I ) * omega[i,])
                       }))      
    }
    samples = data.frame( pi_ij, pi_star, omega )
  }
  
  cn1 = paste0("pi[",outer(1:nrow(data), 1:ncol(data),paste,sep=","),"]")
  if(fixedMargin == "rows"){
    cn2 = paste0("pi[",1:nrow(data),",*]") 
  }else{
    cn2 = paste0("pi[*,",1:ncol(data),"]") 
  }
  cn3 = paste0("omega[",outer(1:nrow(data), 1:ncol(data),paste,sep=","),"]")
  colnames(samples) = c(cn1,cn2,cn3)
  
  return(mcmc(samples))
  
}

sampleHypergeomContingencyNull <- function(prior, data, iterations, noSample = FALSE, ...)
{
    
  if(noSample){
    samples = data.frame(matrix(as.numeric(NA), 1))
  }else{
    samples = data.frame(matrix(0, iterations, 1))
  }
  
  colnames(samples) = c("log.odds.ratio")
  
  return(mcmc(samples))
}


sampleHypergeomContingencyAlt <- function(prior, data, iterations, noSample = FALSE, ...)
{
  if(noSample){
    samples = data.frame(matrix(as.numeric(NA), 1))
  }else{
    stop("Sampling for this model not yet implemented.")
  }
  
  colnames(samples) = c("log.odds.ratio")
  
  return(mcmc(samples))
  
}


