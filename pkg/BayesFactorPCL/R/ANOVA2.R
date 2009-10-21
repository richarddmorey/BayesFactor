library(coda)
library(mvtnorm)
library(MCMCpack)
library(sciplot)

N = 120   # Total Observations
p=5      # Number of nonmean parameters
K = 2    # number of parameter groups

# design matrix
X = cbind(1,rep(c(1,0),each=60),rep(c(0,1),each=60))
X = cbind(X,rep(c(1,0,0),2*20),rep(c(0,1,0),2*20),rep(c(0,0,1),2*20))
                                        #X = cbind(1,rep(c(1,0,0),each=N/p),rep(c(0,1,0),each=N/p),rep(c(0,0,1),each=N/p))

# Number of parameters in each group
Nks = c(2,3)
grps=c(0,inverse.rle(list(values=1:K,lengths=Nks)))

# True parameters
true.beta=matrix(c(3,-.4,.4,-.5,0,.5),ncol=1)

true.sig2 = 1


# Draw data
y = rnorm(N,X%*%true.beta,sqrt(true.sig2))


# reserve space for chains, and init starting values
start.g = rep(1,K)

chain.beta = matrix(NA,nrow=length(true.beta),ncol=M)
chain.g    = matrix(NA,nrow=K,ncol=M)
chain.sig2 = 1:M * NA

chain.beta[,1] = true.beta
chain.sig2[1] = true.sig2
chain.g[,1] = start.g

B = diag(c(0,inverse.rle(list(values=1/start.g,lengths=Nks))))
XtX = tcrossprod(t(X))

dens = matrix(NA,K,M-1)


# length of chain
M=1000


for(m in 2:M)
{
	# Draw g
	for(k in 1:K)
	{
		sub.bet = chain.beta[grps==k,m-1]
		chain.g[k,m] = rinvgamma(1, Nks[k]/2 + .5,.5 + .5*sum(sub.bet^2)/chain.sig2[m-1])
	}


        # draw sigma2
	B = diag(c(0,inverse.rle(list(values=1/chain.g[,m],lengths=Nks))))
        chain.sig2[m] = rinvgamma(1,N/2+p/2,.5*sum((y-X%*%chain.beta[,m-1])^2)+.5*sum(chain.beta[,m-1]^2*diag(B)))


        # draw beta
	Sigma = chain.sig2[m] * solve(XtX + B)
	Mu    = Sigma%*%t(X)%*%y/chain.sig2[m]
	
	chain.beta[,m] = mvrnorm(1,Mu,Sigma)

        # Theorem for conditional distributions
        # of subvectors of multivariate normal vectors
        # See Moser (or other linear models text)
        for(k in 1:K)
	{
          subset = grps==k
          mu1 = Mu[subset]/sqrt(chain.sig2[m])
          mu2 = Mu[!subset]/sqrt(chain.sig2[m])
          known.beta = chain.beta[!subset,m]
          S11 = Sigma[subset,subset]/chain.sig2[m]
          S22 = Sigma[!subset,!subset]/chain.sig2[m]
          S12 = Sigma[subset,!subset]/chain.sig2[m]
          condMean = mu1 + S12%*%solve(S22)%*%(known.beta-mu2)
          condVariance = S11 - S12%*%solve(S22)%*%t(S12)
          dens[k,m-1] = dmvnorm(mu1*0,condMean,condVariance)
          
        }
}


MCMCN = apply(dens,1,effectiveSize)
vars  = apply(dens,1,var)

prior=(1/pi)^Nks

matplot(t(chain.beta))
rowMeans(chain.beta)
tapply(y,rep(1:p,each=N/p),mean)

plot(chain.sig2,ty='l')


BF = rowMeans(dens)/prior
CI = rbind(BF + 1.96*sqrt(vars/MCMCN),BF - 1.96*sqrt(vars/MCMCN))

  
a1 = as.factor(rep(c(1,2),each=60))
a2 = as.factor(rep(c(1,2,3),40))
summary(aov(y~a1+a2))
bargraph.CI(a1,y,a2)

