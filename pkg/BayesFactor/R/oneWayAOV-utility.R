
marginal.g.oneWay = Vectorize(function(g,F,N,J,rscale,log=FALSE, log.const=0)
{
  dfs = (J-1)/(N*J-J)
  omega = (1+(N*g/(dfs*F+1)))/(N*g+1)
  m = log(rscale) - 0.5*log(2*pi) - 1.5*log(g) - rscale^2/(2*g) - (J-1)/2*log(N*g+1) - (N*J-1)/2*log(omega) - log.const
  ifelse(log,m,exp(m))
},"g")

