## ----echo=FALSE,message=FALSE,results='hide'-----------------------------


## ----echo=FALSE,message=FALSE,results='hide'-----------------------------
options(markdown.HTML.stylesheet = 'extra/manual.css')
library(knitr)
library(BayesFactor)
options(BFprogress = FALSE)
bfversion = BFInfo()
session = sessionInfo()[[1]]
rversion = paste(session$version.string," on ",session$platform,sep="")
set.seed(2)

## ------------------------------------------------------------------------
# Create data
x <- rnorm(20)
x[1:10] = x[1:10] + .2
grp = factor(rep(1:2,each=10))

dat = data.frame(x=x,grp=grp)

t.test(x ~ grp, data=dat)

## ------------------------------------------------------------------------
as.vector(ttestBF(formula = x ~ grp, data=dat))
as.vector(anovaBF(x~grp, data=dat))
as.vector(generalTestBF(x~grp, data=dat))

## ------------------------------------------------------------------------
# create some data
id = rnorm(10)
eff = c(-1,1)*1
effCross = outer(id,eff,'+')+rnorm(length(id)*2)
dat = data.frame(x=as.vector(effCross),id=factor(1:10), grp=factor(rep(1:2,each=length(id))))
dat$forReg = as.numeric(dat$grp)-1.5
idOnly = lmBF(x~id, data=dat, whichRandom="id")

summary(aov(x~grp+Error(id/grp),data=dat))


## ------------------------------------------------------------------------
as.vector(lmBF(x ~ grp+id, data=dat, whichRandom="id")/idOnly)
as.vector(lmBF(x ~ forReg+id, data=dat, whichRandom="id")/idOnly)

## ------------------------------------------------------------------------
# create some data
tstat = 3
NTwoSample = 500
effSampleSize = (NTwoSample^2)/(2*NTwoSample)
effSize = tstat/sqrt(effSampleSize)

# One sample
x0 = rnorm(effSampleSize)
x0 = (x0 - mean(x0))/sd(x0) + effSize

t.test(x0)

# Two sample
x1 = rnorm(NTwoSample)
x1 = (x1 - mean(x1))/sd(x1)
x2 = x1 + effSize

t.test(x2,x1)


## ------------------------------------------------------------------------
log(as.vector(ttestBF(x0)))
log(as.vector(ttestBF(x=x1,y=x2)))

## ------------------------------------------------------------------------
# using the data previously defined
t.test(x~grp,data=dat,paired=TRUE)

as.vector(lmBF(x ~ grp+id, data=dat, whichRandom="id")/idOnly)
as.vector(ttestBF(x=dat$x[dat$grp==1],y=dat$x[dat$grp==2],paired=TRUE))

