## ----echo=FALSE,message=FALSE,results='hide'-----------------------------
library(BayesFactor)

options(BFprogress = FALSE)
bfversion = BFInfo()
session = sessionInfo()[[1]]
rversion = paste(session$version.string," on ",session$platform,sep="")


options(markdown.HTML.stylesheet = 'extra/manual.css')
library(knitr)
options(digits=3)
require(graphics)
set.seed(2)

## ----message=FALSE,warning=FALSE-----------------------------------------
library(arm)
library(lme4)


## ------------------------------------------------------------------------
# Number of participants
N <- 20
sig2 <- 1
sig2ID <- 1

# 3x3x3 design, with participant as random factor
effects <- expand.grid(A = c("A1","A2","A3"),
                       B = c("B1","B2","B3"),
                       C = c("C1","C2","C3"),
                       ID = paste("Sub",1:N,sep="")
)
Xdata <- model.matrix(~ A*B*C + ID, data=effects)
beta <- matrix(c(50,
          -.2,.2,
          0,0,
          .1,-.1,
          rnorm(N-1,0,sqrt(sig2ID)),
          0,0,0,0,
          -.1,.1,.1,-.1,
          0,0,0,0,
          0,0,0,0,0,0,0,0),
               ncol=1)
effects$y = rnorm(Xdata%*%beta,Xdata%*%beta,sqrt(sig2))

## ------------------------------------------------------------------------
# Typical repeated measures ANOVA
summary(fullaov <- aov(y ~ A*B*C + Error(ID/(A*B*C)),data=effects))

## ----fig.width=10,fig.height=4-------------------------------------------
mns <- tapply(effects$y,list(effects$A,effects$B,effects$C),mean)
stderr = sqrt((sum(resid(fullaov[[3]])^2)/fullaov[[3]]$df.resid)/N)

par(mfrow=c(1,3),cex=1.1)
for(i in 1:3){
  matplot(mns[,,i],xaxt='n',typ='b',xlab="A",main=paste("C",i), 
          ylim=range(mns)+c(-1,1)*stderr,ylab="y")
  axis(1,at=1:3,lab=1:3)
  segments(1:3 + mns[,,i]*0,mns[,,i] + stderr,1:3 + mns[,,i]*0,mns[,,i] - stderr,col=rgb(0,0,0,.3))
}

## ------------------------------------------------------------------------
t.is = system.time(bfs.is <- anovaBF(y ~ A*B*C + ID, data = effects, 
                                     whichRandom="ID")
)
t.la = system.time(bfs.la <- anovaBF(y ~ A*B*C + ID, data = effects, 
                                     whichRandom="ID",
                                     method = "laplace")
)

## ----fig.width=6,fig.height=6--------------------------------------------
t.is
t.la

plot(log(extractBF(sort(bfs.is))$bf),log(extractBF(sort(bfs.la))$bf),
     xlab="Default Sampler",ylab="Laplace approximation",
     pch=21,bg=rgb(0,0,1,.2),col="black",asp=TRUE,cex=1.2)
abline(0,1)

bfs.is

## ----message=FALSE-------------------------------------------------------
chains <- lmBF(y ~ A + B + C + ID, data=effects, whichRandom = "ID", posterior=TRUE, iterations=10000)

lmerObj <- lmer(y ~ A + B + C + (1|ID), data=effects)
# Use arm function sim() to sample from posterior
chainsLmer = sim(lmerObj,n.sims=10000)

## ------------------------------------------------------------------------
BF.sig2 <- chains[,colnames(chains)=="sig2"]
AG.sig2 <- (chainsLmer@sigma)^2
qqplot(log(BF.sig2),log(AG.sig2),pch=21,bg=rgb(0,0,1,.2),
       col=NULL,asp=TRUE,cex=1,xlab="BayesFactor samples",
       ylab="arm samples",main="Posterior samples of\nerror variance")
abline(0,1)

## ------------------------------------------------------------------------
AG.raneff <- chainsLmer@ranef$ID[,,1]
BF.raneff <-  chains[,grep('ID-',colnames(chains),fixed='TRUE')]
plot(colMeans(BF.raneff),colMeans(AG.raneff),pch=21,bg=rgb(0,0,1,.2),col="black",asp=TRUE,cex=1.2,xlab="BayesFactor estimate",ylab="arm estimate",main="Random effect posterior means")
abline(0,1)


## ----tidy=FALSE----------------------------------------------------------
AG.fixeff <- chainsLmer@fixef
BF.fixeff <-  chains[,1:10]

# Adjust AG results from reference cell to sum to 0
Z = c(1,  1/3,  1/3,  1/3,  1/3,  1/3,  1/3,
      0, -1/3, -1/3,    0,    0,    0,    0,
      0,  2/3, -1/3,    0,    0,    0,    0,
      0, -1/3,  2/3,    0,    0,    0,    0,
      0,     0,   0, -1/3, -1/3,    0,    0,
      0,     0,   0,  2/3, -1/3,    0,    0,
      0,     0,   0, -1/3,  2/3,    0,    0,
      0,     0,   0,    0,    0, -1/3, -1/3,
      0,     0,   0,    0,    0,  2/3, -1/3,
      0,     0,   0,    0,    0, -1/3,  2/3)
dim(Z) = c(7,10)
Z = t(Z)

AG.fixeff2 = t(Z%*%t(AG.fixeff))

## Our grand mean has heavier tails
qqplot(BF.fixeff[,1],AG.fixeff2[,1],pch=21,bg=rgb(0,0,1,.2),col=NULL,asp=TRUE,cex=1,xlab="BayesFactor estimate",ylab="arm estimate",main="Grand mean posterior samples")
abline(0,1)

plot(colMeans(BF.fixeff[,-1]),colMeans(AG.fixeff2[,-1]),pch=21,bg=rgb(0,0,1,.2),col="black",asp=TRUE,cex=1.2,xlab="BayesFactor estimate",ylab="arm estimate",main="Fixed effect posterior means")
abline(0,1)

## Compare posterior standard deviations
BFsd = apply(BF.fixeff[,-1],2,sd)
AGsd = apply(AG.fixeff2[,-1],2,sd)
plot(sort(AGsd/BFsd),pch=21,bg=rgb(0,0,1,.2),col="black",cex=1.2,ylab="Ratio of posterior standard deviations (arm/BF)",xlab="Fixed effect index")
## AG estimates are slightly larger, consistent with sig2 estimates
## probably due to prior


## ----message=FALSE,warning=FALSE-----------------------------------------
library(languageR)
library(xtable)

## ------------------------------------------------------------------------
data(primingHeidPrevRT)

primingHeidPrevRT$lRTmin1 <- log(primingHeidPrevRT$RTmin1)

###Frequentist 

lr4 <- lmer(RT ~ Condition + (1|Word)+ (1|Subject) + lRTmin1 + RTtoPrime + ResponseToPrime + ResponseToPrime*RTtoPrime +BaseFrequency ,primingHeidPrevRT)
# Get rid rid of some outlying response times
INDOL <- which(scale(resid(lr4)) < 2.5)
primHeidOL <- primingHeidPrevRT[INDOL,]

## ------------------------------------------------------------------------
# Center continuous variables
primHeidOL$BaseFrequency <- primHeidOL$BaseFrequency - mean(primHeidOL$BaseFrequency)
primHeidOL$lRTmin1 <- primHeidOL$lRTmin1 - mean(primHeidOL$lRTmin1)
primHeidOL$RTtoPrime <- primHeidOL$RTtoPrime - mean(primHeidOL$RTtoPrime)

## ------------------------------------------------------------------------
# LMER
lr4b <- lmer(  RT ~ Condition + ResponseToPrime +  (1|Word)+ (1|Subject) + lRTmin1 + RTtoPrime + ResponseToPrime*RTtoPrime + BaseFrequency , primHeidOL)
# BayesFactor
B5out <- lmBF( RT ~ Condition + ResponseToPrime +     Word +    Subject  + lRTmin1 + RTtoPrime + ResponseToPrime*RTtoPrime + BaseFrequency  , primHeidOL , whichRandom = c("Word", "Subject"),  posterior = TRUE, iteration = 50000,columnFilter=c("Word","Subject"))

lmerEff <- fixef(lr4b)
bfEff <- colMeans(B5out[,1:10])

## ----results='asis'------------------------------------------------------
print(xtable(cbind("lmer fixed effects"=names(lmerEff))), type='html')

## ----tidy=FALSE----------------------------------------------------------
# Adjust lmer results from reference cell to sum to 0
Z = c(1,   1/2, 1/2,    0,    0,    0,    0,
      0,  -1/2,   0,    0,    0,    0,    0,
      0,   1/2,   0,    0,    0,    0,    0,
      0,     0,-1/2,    0,    0,    0,    0,
      0,     0, 1/2,    0,    0,    0,    0,
      0,     0,   0,    1,    0,    0,    0,
      0,     0,   0,    0,    1,    0,  1/2,
      0,     0,   0,    0,    0,    1,    0,
      0,     0,   0,    0,    0,    0, -1/2,
      0,     0,   0,    0,    0,    0,  1/2)
dim(Z) = c(7,10)
Z = t(Z)

# Do reparameterization by pre-multimplying the parameter vector by Z
reparLmer <- Z %*% matrix(lmerEff,ncol=1)

# put results in data.frame for comparison
sideBySide <- data.frame(BayesFactor=bfEff,lmer=reparLmer)

## ----results='asis'------------------------------------------------------
print(xtable(sideBySide,digits=4), type='html')

## ------------------------------------------------------------------------
# Notice Bayesian shrinkage
par(cex=1.5)
plot(sideBySide[-1,],pch=21,bg=rgb(0,0,1,.2),col="black",asp=TRUE,cex=1.2, main="fixed effects\n (excluding grand mean)")
abline(0,1, lty=2)

