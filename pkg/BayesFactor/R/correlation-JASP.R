## The code in this file is from the JASP project
## (https://github.com/jasp-stats/jasp-desktop/blob/development/JASP-Engine/JASP/R/correlationbayesian.R)
## and written by Alexander Ly (Alexander.Ly.NL@gmail.com)

.bf10Exact <- function(n, r, kappa=1) {
  # Ly et al 2015
  # This is the exact result with symmetric beta prior on rho
  # with parameter alpha. If kappa = 1 then uniform prior on rho
  #
  #
  if (n <= 2){
    return(list(bf = 0, properror=0,method="exact"))
  } else if (any(is.na(r))){
    return(list(bf = NA, properror=NA,method=NA))

  }
  # TODO: use which
  check.r <- abs(r) >= 1 # check whether |r| >= 1
  if (kappa >= 1 && n > 2 && check.r) {
    return(list(bf = Inf, properror=0,method="exact"))

  }
  #log.hyper.term <- log(hypergeo::hypergeo(((n-1)/2), ((n-1)/2), ((n+2/kappa)/2), r^2))
  log.hyper.term <- Re(
    genhypergeo_series_pos(U=c((n-1)/2, (n-1)/2),
                           L=((n+2/kappa)/2), z=r^2,
                           0, 2000, TRUE, TRUE, FALSE)
    )
  log.result <- (1-2/kappa)*log(2)+0.5*log(pi)-lbeta(1/kappa, 1/kappa)+
    lgamma((n+2/kappa-1)/2)-lgamma((n+2/kappa)/2)+log.hyper.term
  real.result <- log.result
  return(list(bf = real.result, properror=0,method="exact"))
}

# 2.2 Two-sided secondary Bayes factor
.bf10JeffreysIntegrate <- function(n, r, kappa=1) {
  # Jeffreys' test for whether a correlation is zero or not
  # Jeffreys (1961), pp. 289-292
  # This is the exact result, see EJ
  ##
  if (n <= 2){
    return(list(bf = 0, properror=0,method="exact"))
  } else if ( any(is.na(r)) ){
    return(list(bf = NA, properror=NA,method=NA))

  }

  # TODO: use which
  if (n > 2 && abs(r)==1) {
    return(list(bf = Inf, properror=0,method="exact"))

  }
  hyper.term <- Re(
    genhypergeo_series_pos(U=c((2*n-3)/4, (2*n-1)/4),
                           L=(n+2/kappa)/2, z=r^2,
                           0, 2000, TRUE, TRUE, FALSE)
    )
  log.term <- lgamma((n+2/kappa-1)/2)-lgamma((n+2/kappa)/2)-lbeta(1/kappa, 1/kappa)
  result <- .5*log(pi) + (1-2/kappa)*log(2) + log.term + hyper.term
  return(list(bf = result, properror=NA, method="Jeffreys' approximation"))
}

# 2.3 Two-sided third Bayes factor
.bfCorNumerical <- function(n, r, kappa=1, lowerRho=-1, upperRho=1, approx = TRUE) {
  # Numerically integrate Jeffreys approximation of the likelihood

  if(approx){
    likeFun = jeffreys_approx_corr
  }else{
    likeFun = function(rho,n,r){
      hFunc(rho, n, r, FALSE, 2000)
    }
  }

  log.const = likeFun(r,n,r)

  integrand <- Vectorize(function(rho){
    exp(likeFun(rho, n, r) +
      dbeta(rho/2+.5,1/kappa,1/kappa,log=TRUE) - log(2) - log.const)
  },"rho")
  some.integral <- try(integrate(integrand, lowerRho, upperRho))

  if (is(some.integral, "try-error")) {
    return(NULL)
  }

  if (some.integral$message=="OK"){
    some.integral$value = exp(log(some.integral$value) + log.const)
    return(some.integral)
  } else {
    return(NULL)
  }
}

.bf10Numerical <- function(n, r, kappa=1, lowerRho=-1, upperRho=1,approx=TRUE) {
  # Jeffreys' test for whether a correlation is zero or not
  # Jeffreys (1961), pp. 289-292
  # This is a numerical approximation for .bf10JeffreysIntegrate,
  # when it explodes
  # #
  # TODO: 1. check for n=1, n=2, as r is then undefined
  # 2. check for r=1, r=-1
  #
  # TODO: REMOVE ALL NUMERICAL STUFF
  if ( any(is.na(r)) ){
    return(list(bf = NA, properror=NA, method=NA))
  }
  # TODO: use which
  if (n > 2 && abs(r)==1) {
    return(list(bf = Inf, properror=0, method="exact"))
  }

  # TODO: be very careful here, might integrate over non-finite function
  jeffreysNumericalIntegrate <- .bfCorNumerical(n, r, kappa, lowerRho=-1, upperRho=1,approx)

  if (is.null(jeffreysNumericalIntegrate)){
    return(list(bf = NA, properror=NA, method=NA))
  }else if(jeffreysNumericalIntegrate$value < 0){
    return(list(bf = NA, properror=NA, method=NA))
  } else if (jeffreysNumericalIntegrate$value >= 0){
    # jeffreys numerical integrate success
    log.bf = log(jeffreysNumericalIntegrate$value)
    err = jeffreysNumericalIntegrate$abs.error
    prop.error = exp(log(err) - log.bf)
    return(list(bf = log.bf, properror=err, method="quadrature"))
  } else {
    # NO IDEA, EVERYTHING FAILED :(
    return(list(bf = NA, properror=NA, method=NA))
  }
  return(list(bf = NA, properror=NA, method=NA))
}

