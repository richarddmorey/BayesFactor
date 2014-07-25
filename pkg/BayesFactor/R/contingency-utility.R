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
  
  lbf<-lgamma(sum(a)-(I-1)*(J-1))-((I-1)*(J-1)*log(1+1/b))-
    lgamma(sum(y)+sum(a)-(I-1)*(J-1))-ldirich(yr + ar-(J- 1)*or) +
    ldirich( ar - (J - 1) * or)-ldirich(yc + ac - (I - 1) * oc)+
    ldirich( ac - (I - 1) * oc)- ldirich1(y,a)
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
    ldirich(ac - (I - 1) * oc) - ldirich(c(a)) - ldirich(yr + 
                                                           ar - (J - 1) * or) - ldirich(yc + ac - (I - 1) * oc)
    
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
  lbf <- ldirich(ac - (I - 1) * oc) + ldirich(ar) + ldirich(c(y) + c(a))-
    ldirich(yc + ac - (I - 1) * oc)- ldirich(yr + ar) - ldirich(c(a))
  
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
  lbf<-ldirich(c(y)+c(a))+pcomb(yc)+pcomb(yr)-ldirich(c(a))-sumg(y,a)
  
  return(lbf)  
}

###########################
## Sampling
###########################

sampleContingency <- function(model, type, prior, data = data, iterations = iterations, ...)                  
{
  stop("Sampling for contingency tables not yet implemented.")
}


