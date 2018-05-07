library(hypergeo)

A <- function(t, n, nu, mu.delta, g) {
  
  Re(hypergeo::genhypergeo(U = (nu + 1)/2, L = 1/2,
                           z = mu.delta^2*t^2/
                             (2*(1/n + g)*((1 + n*g)*nu + t^2))))
  
}

B <- function(t, n, nu, mu.delta, g) {
  
  out <- mu.delta*t/sqrt(1/2*(1/n + g)*((1 + n*g)*nu + t^2)) *
    exp(lgamma((nu + 2)/2) - lgamma((nu + 1)/2)) *
    Re(hypergeo::genhypergeo(U = (nu + 2)/2, L = 3/2,
                             z = mu.delta^2*t^2/
                               (2*(1/n + g)*((1 + n*g)*nu + t^2))))
  
  return(out)
  
}


C <- function(delta, t, n, nu) {
  
  Re(hypergeo::genhypergeo(U = (nu + 1)/2, L = 1/2,
                           z = n*t^2*delta^2/(2*(nu + t^2))))
  
}

D <- function(delta, t, n, nu) {
  
  out <- t*delta*sqrt(2*n/(nu + t^2))*
    exp(lgamma((nu + 2)/2) - lgamma((nu + 1)/2))*
    Re(hypergeo::genhypergeo(U = (nu + 2)/2, L = 3/2,
                             z = n*t^2*delta^2/(2*(nu + t^2))))
  
  return(out)
  
}

term_normalprior <- function(t, n, nu, mu.delta, g) {
  
  (1 + n*g)^(-1/2) * exp(-mu.delta^2/(2*(1/n + g))) *
    (1 + t^2/(nu*(1 + n*g)))^(-(nu + 1)/2) *
    (A(t, n, nu, mu.delta, g) + B(t, n, nu, mu.delta, g))
  
}

integrand <- function(g, t, n, nu, mu.delta, r, kappa) {
  
  tmp <- term_normalprior(t = t, n = n, nu = nu, mu.delta = mu.delta, g = g)
  pg_log <- kappa/2*(2*log(r) + log(kappa/2)) - lgamma(kappa/2) -
    (kappa/2 + 1)*log(g) - r^2*kappa/(2*g)
  pg <- exp(pg_log)
  out <- tmp*pg
  
  return(out) 
  
}

dtss <- function(delta, mu.delta, r, kappa, log = FALSE) {
  
  out <- - log(r) + lgamma((kappa + 1)/2) - .5*(log(pi) + log(kappa)) - 
    lgamma(kappa/2) - (kappa + 1)/2 * log(1 + ((delta - mu.delta)/r)^2/kappa)
  
  if ( ! log)
    out <- exp(out)
  
  return(out)
  
}

posterior_t_tmp <- function(delta, t, ny, nx = NULL, independentSamples = FALSE,
                            prior.location, prior.scale, prior.df,
                            rel.tol = .Machine$double.eps^0.25) {
  
  neff <- ifelse(independentSamples, ny*nx/(ny + nx), ny)
  nu <- ifelse(independentSamples, ny + nx - 2, ny - 1)
  
  mu.delta <- prior.location
  r <- prior.scale
  kappa <- prior.df
  
  numerator <- exp(-neff/2*delta^2)*(1 + t^2/nu)^(-(nu + 1)/2)*
    (C(delta, t, neff, nu) + D(delta, t, neff, nu))*
    dtss(delta, mu.delta, r, kappa)
  
  denominator <- integrate(integrand, lower = 0, upper = Inf,
                           t = t, n = neff, nu = nu, mu.delta = mu.delta,
                           r = r, kappa = kappa, rel.tol = rel.tol)$value
  
  out <- numerator/denominator
  
  if ( is.na(out))
    out <- 0
  
  return(out)
  
}

posterior_t <- Vectorize(posterior_t_tmp, "delta")

cdf_t <- function(x, t, ny, nx = NULL, independentSamples = FALSE,
                  prior.location, prior.scale, prior.df) {
  
  integrate(posterior_t, lower = -Inf, upper = x, t = t, ny = ny, nx = nx,
            independentSamples = independentSamples,
            prior.location = prior.location, prior.scale = prior.scale,
            prior.df = prior.df)$value
  
}

quantile_t <- function(q, t, ny, nx = NULL,
                       independentSamples = FALSE,
                       prior.location, prior.scale,
                       prior.df, tol = 0.0001, max.iter = 100) {
  
  # compute quantiles via Newton-Raphson method
  
  x.cur <- Inf
  # get reasonable starting value
  delta <- seq(-2, 2, length.out = 400)
  dens <- posterior_t(delta, t = t, ny = ny, nx = nx,
                      independentSamples = independentSamples,
                      prior.location = prior.location,
                      prior.scale = prior.scale,
                      prior.df = prior.df)
  x.new <- delta[which.max(dens)]
  i <- 1
  
  while (abs(x.cur - x.new) > tol && i < max.iter) {
    
    x.cur <- x.new
    x.new <- x.cur - (cdf_t(x.cur, t = t, ny = ny, nx = nx,
                            independentSamples = independentSamples,
                            prior.location = prior.location,
                            prior.scale = prior.scale,
                            prior.df = prior.df) - q)/
      posterior_t(x.cur, t = t, ny = ny, nx = nx,
                  independentSamples = independentSamples,
                  prior.location = prior.location, prior.scale = prior.scale,
                  prior.df = prior.df)
    i <- i + 1
    
  }
  
  return(x.new)
  
}

ciPlusMedian_t <- function(t, ny, nx = NULL, independentSamples = FALSE,
                           prior.location, prior.scale, prior.df,
                           ci = .95, type = "two-sided", tol = 0.0001, max.iter = 100) {
  
  lower <- (1 - ci)/2
  upper <- ci + (1 - ci)/2
  med <- .5
  
  postAreaSmaller0 <- cdf_t(x = 0, t = t, ny = ny, nx = nx,
                            independentSamples = independentSamples,
                            prior.location = prior.location,
                            prior.scale = prior.scale, prior.df = prior.df)
  
  if (type == "plus-sided") {
    
    lower <- postAreaSmaller0 + (1 - postAreaSmaller0)*lower
    upper <- postAreaSmaller0 + (1 - postAreaSmaller0)*upper
    med <- postAreaSmaller0 + (1 - postAreaSmaller0)*med
    
  } else if (type == "min-sided") {
    
    lower <- postAreaSmaller0*lower
    upper <- postAreaSmaller0*upper
    med <- postAreaSmaller0*med
    
  }
  
  ciLower <- quantile_t(lower, t = t, ny = ny, nx = nx,
                        independentSamples = independentSamples,
                        prior.location = prior.location,
                        prior.scale = prior.scale,
                        prior.df = prior.df)
  ciUpper <- quantile_t(upper, t = t, ny = ny, nx = nx,
                        independentSamples = independentSamples,
                        prior.location = prior.location,
                        prior.scale = prior.scale,
                        prior.df = prior.df)
  median <- quantile_t(med, t = t, ny = ny, nx = nx,
                       independentSamples = independentSamples,
                       prior.location = prior.location,
                       prior.scale = prior.scale,
                       prior.df = prior.df)
  
  return(list(ciLower = ciLower, median = median, ciUpper = ciUpper))
  
}

posterior_normal_tmp <- function(delta, t, ny, nx = NULL,
                                 independentSamples = FALSE, prior.mean,
                                 prior.variance,
                                 rel.tol = .Machine$double.eps^0.25) {
  
  neff <- ifelse(independentSamples, ny*nx/(ny + nx), ny)
  nu <- ifelse(independentSamples, ny + nx - 2, ny - 1)
  
  mu.delta <- prior.mean
  g <- prior.variance
  
  numerator <- exp(-neff/2*delta^2)*(1 + t^2/nu)^(-(nu + 1)/2)*
    (C(delta, t, neff, nu) + D(delta, t, neff, nu))*
    dnorm(delta, mu.delta, sqrt(g))
  
  denominator <- term_normalprior(t = t, n = neff, nu = nu,
                                  mu.delta = mu.delta, g = g)
  
  out <- numerator/denominator
  
  if ( is.na(out))
    out <- 0
  
  return(out)
  
}

posterior_normal <- Vectorize(posterior_normal_tmp, "delta")

cdf_normal <- function(x, t, ny, nx = NULL, independentSamples = FALSE,
                       prior.mean, prior.variance) {
  
  integrate(posterior_normal, lower = -Inf, upper = x, t = t, ny = ny, nx = nx,
            independentSamples = independentSamples,
            prior.mean = prior.mean, prior.variance = prior.variance)$value
  
}

quantile_normal <- function(q, t, ny, nx = NULL,
                            independentSamples = FALSE,
                            prior.mean, prior.variance,
                            tol = 0.0001, max.iter = 100) {
  
  # compute quantiles via Newton-Raphson method
  
  x.cur <- Inf
  # get reasonable start value
  delta <- seq(-2, 2, length.out = 400)
  dens <- posterior_normal(delta, t = t, ny = ny, nx = nx,
                           independentSamples = independentSamples,
                           prior.mean = prior.mean,
                           prior.variance = prior.variance)
  x.new <- delta[which.max(dens)]
  i <- 1
  
  while (abs(x.cur - x.new) > tol && i < max.iter) {
    
    x.cur <- x.new
    x.new <- x.cur - (cdf_normal(x.cur, t = t, ny = ny, nx = nx,
                                 independentSamples = independentSamples,
                                 prior.mean = prior.mean,
                                 prior.variance = prior.variance) - q)/
      posterior_normal(x.cur, t = t, ny = ny, nx = nx,
                       independentSamples = independentSamples,
                       prior.mean = prior.mean, prior.variance = prior.variance)
    i <- i + 1
    
  }
  
  return(x.new)
  
}

ciPlusMedian_normal <- function(t, ny, nx = NULL, independentSamples = FALSE,
                                prior.mean, prior.variance, ci = .95,
                                type = "two-sided", tol = 0.0001, max.iter = 100) {
  
  lower <- (1 - ci)/2
  upper <- ci + (1 - ci)/2
  med <- .5
  
  postAreaSmaller0 <- cdf_normal(x = 0, t = t, ny = ny, nx = nx,
                                 independentSamples = independentSamples,
                                 prior.mean = prior.mean,
                                 prior.variance = prior.variance)
  
  if (type == "plus-sided") {
    
    lower <- postAreaSmaller0 + (1 - postAreaSmaller0)*lower
    upper <- postAreaSmaller0 + (1 - postAreaSmaller0)*upper
    med <- postAreaSmaller0 + (1 - postAreaSmaller0)*med
    
  } else if (type == "min-sided") {
    
    lower <- postAreaSmaller0*lower
    upper <- postAreaSmaller0*upper
    med <- postAreaSmaller0*med
    
  }
  
  ciLower <- quantile_normal(lower, t = t, ny = ny, nx = nx,
                             independentSamples = independentSamples,
                             prior.mean = prior.mean,
                             prior.variance = prior.variance)
  ciUpper <- quantile_normal(upper, t = t, ny = ny, nx = nx,
                             independentSamples = independentSamples,
                             prior.mean = prior.mean,
                             prior.variance = prior.variance)
  median <- quantile_normal(med, t = t, ny = ny, nx = nx,
                            independentSamples = independentSamples,
                            prior.mean = prior.mean,
                            prior.variance = prior.variance)
  
  return(list(ciLower = ciLower, median = median, ciUpper = ciUpper))
  
}

bf10_t <- function(t, ny, nx = NULL, independentSamples = FALSE, prior.location,
                   prior.scale, prior.df, rel.tol = .Machine$double.eps^0.25) {
  
  neff <- ifelse(independentSamples, ny*nx/(ny + nx), ny)
  nu <- ifelse(independentSamples, ny + nx - 2, ny - 1)
  
  mu.delta <- prior.location
  r <- prior.scale
  kappa <- prior.df
  numerator <- integrate(integrand, lower = 0, upper = Inf,
                         t = t, n = neff, nu = nu, mu.delta = mu.delta,
                         r = r, kappa = kappa,
                         rel.tol = rel.tol)$value
  denominator <- (1 + t^2/nu)^(-(nu + 1)/2)
  
  BF10 <- numerator/denominator
  priorAreaSmaller0 <- integrate(dtss, lower = -Inf, upper = 0,
                                 mu.delta = prior.location, r = prior.scale,
                                 kappa = prior.df)$value
  postAreaSmaller0 <- cdf_t(x = 0, t = t, ny = ny, nx = nx,
                            independentSamples = independentSamples,
                            prior.location = prior.location,
                            prior.scale = prior.scale, prior.df = prior.df)
  BFmin1 <- postAreaSmaller0/priorAreaSmaller0
  BFplus1 <- (1 - postAreaSmaller0)/(1 - priorAreaSmaller0)
  BFmin0 <- BFmin1 * BF10
  BFplus0 <- BFplus1 * BF10
  
  return(list(BF10 = BF10, BFplus0 = BFplus0, BFmin0 = BFmin0))
  
}

bf10_normal <- function(t, ny, nx = NULL, independentSamples = FALSE,
                        prior.mean, prior.variance) {
  
  neff <- ifelse(independentSamples, ny*nx/(ny + nx), ny)
  nu <- ifelse(independentSamples, ny + nx - 2, ny - 1)
  
  mu.delta <- prior.mean
  g <- prior.variance
  numerator <- term_normalprior(t = t, n = neff, nu  = nu,
                                mu.delta = mu.delta, g = g)
  denominator <- (1 + t^2/nu)^(-(nu + 1)/2)
  
  BF10 <- numerator/denominator
  priorAreaSmaller0 <- pnorm(0, mean = prior.mean, sd = sqrt(prior.variance))
  postAreaSmaller0 <- cdf_normal(x = 0, t = t, ny = ny, nx = nx,
                                 independentSamples = independentSamples,
                                 prior.mean = prior.mean,
                                 prior.variance = prior.variance)
  BFmin1 <- postAreaSmaller0/priorAreaSmaller0
  BFplus1 <- (1 - postAreaSmaller0)/(1 - priorAreaSmaller0)
  BFmin0 <- BFmin1 * BF10
  BFplus0 <- BFplus1 * BF10
  
  return(list(BF10 = BF10, BFplus0 = BFplus0, BFmin0 = BFmin0))
  
}
