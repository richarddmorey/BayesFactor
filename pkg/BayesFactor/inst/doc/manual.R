## ----echo=FALSE,message=FALSE,results='hide'-----------------------------
options(markdown.HTML.stylesheet = 'extra/manual.css')
library(knitr)
options(digits=3)
require(graphics)
set.seed(2)

## ----message=FALSE-------------------------------------------------------
library(BayesFactor)

## ----echo=FALSE,message=FALSE,results='hide'-----------------------------
options(BFprogress = FALSE)
bfversion = BFInfo()
session = sessionInfo()[[1]]
rversion = paste(session$version.string," on ",session$platform,sep="")

## ----onesampdata---------------------------------------------------------
data(sleep)

## Compute difference scores
diffScores = sleep$extra[1:10] - sleep$extra[11:20]

## Traditional two-tailed t test
t.test(diffScores)

## ----onesampt------------------------------------------------------------
bf = ttestBF(x = diffScores)
## Equivalently:
## bf = ttestBF(x = sleep$extra[1:10],y=sleep$extra[11:20], paired=TRUE)
bf

## ----recip---------------------------------------------------------------
1 / bf

## ----tsamp---------------------------------------------------------------
chains = posterior(bf, iterations = 1000)
summary(chains)

## ----tsamplplot,fig.width=10---------------------------------------------
chains2 = recompute(chains, iterations = 10000)
plot(chains2[,1:2])

## ----onesamptinterval----------------------------------------------------
bfInterval = ttestBF(x = diffScores, nullInterval=c(-Inf,0))
bfInterval

## ----onesampledivide-----------------------------------------------------
bfInterval[1] / bfInterval[2]

## ----onesampcat----------------------------------------------------------
allbf = c(bf, bfInterval)
allbf

## ----plotonesamp,fig.width=10,fig.height=5-------------------------------
plot(allbf)

## ----onesamplist---------------------------------------------------------
bfmat = allbf / allbf
bfmat

## ----onesamplist2--------------------------------------------------------
bfmat[,2]
bfmat[1,]

## ----onesamplist3--------------------------------------------------------
bfmat[,1:2]
t(bfmat[,1:2])

## ----twosampledata-------------------------------------------------------
data(chickwts)

## Restrict to two groups
chickwts = chickwts[chickwts$feed %in% c("horsebean","linseed"),]
## Drop unused factor levels
chickwts$feed = factor(chickwts$feed)

## Plot data
plot(weight ~  feed, data = chickwts, main = "Chick weights")

## ------------------------------------------------------------------------
## traditional t test
t.test(weight ~ feed, data = chickwts, var.eq=TRUE)

## ----twosamplet----------------------------------------------------------
## Compute Bayes factor
bf = ttestBF(formula = weight ~ feed, data = chickwts)
bf

## ----twosampletsamp,fig.width=10-----------------------------------------
chains = posterior(bf, iterations = 10000)
plot(chains[,1:4])

## ----bemdata-------------------------------------------------------------
## Bem's t statistics from four selected experiments
t = c(-.15, 2.39, 2.42, 2.43)
N = c(100, 150, 97, 99)

## ----bemanalysis1--------------------------------------------------------
bf = meta.ttestBF(t=t, n1=N, nullInterval=c(0,Inf), rscale=1)
bf

## ----bemposterior,fig.width=10-------------------------------------------
## Do analysis again, without nullInterval restriction
bf = meta.ttestBF(t=t, n1=N, rscale=1)

## Obtain posterior samples
chains = posterior(bf, iterations = 10000)
plot(chains)

## ----fixeddata,fig.width=10,fig.height=5---------------------------------
data(ToothGrowth)

## Example plot from ?ToothGrowth

coplot(len ~ dose | supp, data = ToothGrowth, panel = panel.smooth,
       xlab = "ToothGrowth data: length vs dose, given type of supplement")

## Treat dose as a factor
ToothGrowth$dose = factor(ToothGrowth$dose)
levels(ToothGrowth$dose) = c("Low", "Medium", "High")

summary(aov(len ~ supp*dose, data=ToothGrowth))

## ------------------------------------------------------------------------
bf = anovaBF(len ~ supp*dose, data=ToothGrowth)
bf

## ----fixedbf,fig.width=10,fig.height=5-----------------------------------
plot(bf[3:4] / bf[2])

## ------------------------------------------------------------------------
bf = anovaBF(len ~ supp*dose, data=ToothGrowth, whichModels="top")
bf

## ------------------------------------------------------------------------
bfMainEffects = lmBF(len ~ supp + dose, data = ToothGrowth)
bfInteraction = lmBF(len ~ supp + dose + supp:dose, data = ToothGrowth)
## Compare the two models
bf = bfInteraction / bfMainEffects
bf

## ------------------------------------------------------------------------
newbf = recompute(bf, iterations = 500000)
newbf

## ------------------------------------------------------------------------
## Sample from the posterior of the full model
chains = posterior(bfInteraction, iterations = 10000)
## 1:13 are the only "interesting" parameters
summary(chains[,1:13])

## ------------------------------------------------------------------------
plot(chains[,4:6])

## ------------------------------------------------------------------------
data(puzzles)

## ----puzzlesplot,fig.width=7,fig.height=5,echo=FALSE---------------------
## plot the data
aovObj = aov(RT ~ shape*color + Error(ID/(shape*color)), data=puzzles)

matplot(t(matrix(puzzles$RT,12,4)),ty='b',pch=19,lwd=1,lty=1,col=rgb(0,0,0,.2), ylab="Completion time", xlab="Condition",xaxt='n')
axis(1,at=1:4,lab=c("round&mono","square&mono","round&color","square&color"))
mns = tapply(puzzles$RT,list(puzzles$color,puzzles$shape),mean)[c(2,4,1,3)]
points(1:4,mns,pch=22,col="red",bg=rgb(1,0,0,.6),cex=2)
# within-subject standard error, uses MSE from ANOVA
stderr = sqrt(sum(aovObj[[5]]$residuals^2)/11)/sqrt(12)
segments(1:4,mns + stderr,1:4,mns - stderr,col="red")

## ------------------------------------------------------------------------
summary(aov(RT ~ shape*color + Error(ID/(shape*color)), data=puzzles))

## ----tidy=FALSE----------------------------------------------------------
bf = anovaBF(RT ~ shape*color + ID, data = puzzles, 
             whichRandom="ID")

## ------------------------------------------------------------------------
bf

## ----testplot,fig.width=10,fig.height=5----------------------------------
plot(bf)

## ------------------------------------------------------------------------
bfWithoutID = lmBF(RT ~ shape*color, data = puzzles)
bfWithoutID

## ------------------------------------------------------------------------
bfOnlyID = lmBF(RT ~ ID, whichRandom="ID",data = puzzles)
bf2 = bfWithoutID / bfOnlyID
bf2

## ------------------------------------------------------------------------
bfall = c(bf,bf2)

## ------------------------------------------------------------------------
bf[4] / bf2

## ----regressData---------------------------------------------------------
data(attitude)

## Traditional multiple regression analysis
lmObj = lm(rating ~ ., data = attitude)
summary(lmObj)

## ----regressAll----------------------------------------------------------
bf = regressionBF(rating ~ ., data = attitude)
length(bf)

## ----regressSelect-------------------------------------------------------
## Choose a specific model
bf["privileges + learning + raises + critical + advance"]
## Best 6 models
head(bf, n=6)
## Worst 4 models
tail(bf, n=4)

## ----regressSelectwhichmax,eval=FALSE------------------------------------
#  ## which model index is the best?
#  which.max(bf)

## ----regressSelectwhichmaxFake,echo=FALSE--------------------------------
## which model index is the best?
BayesFactor::which.max(bf)

## ----regressSelect2------------------------------------------------------

## Compare the 5 best models to the best
bf2 = head(bf) / max(bf)
bf2
plot(bf2)

## ----regresstop, fig.width=10, fig.height=5------------------------------
bf = regressionBF(rating ~ ., data = attitude, whichModels = "top")
## The seventh model is the most complex
bf
plot(bf)

## ----regressbottom, fig.width=10, fig.height=5---------------------------
bf = regressionBF(rating ~ ., data = attitude, whichModels = "bottom")
plot(bf)

## ----lmregress1----------------------------------------------------------
complaintsOnlyBf = lmBF(rating ~ complaints, data = attitude) 
complaintsLearningBf = lmBF(rating ~ complaints + learning, data = attitude) 
## Compare the two models
complaintsOnlyBf / complaintsLearningBf

## ----lmposterior---------------------------------------------------------
chains = posterior(complaintsLearningBf, iterations = 10000)
summary(chains)

## ----lmregressclassical--------------------------------------------------
summary(lm(rating ~ complaints + learning, data = attitude))

## ----echo=FALSE,results='hide'-------------------------------------------
rm(ToothGrowth)

## ----GLMdata-------------------------------------------------------------
data(ToothGrowth)

# model log2 of dose instead of dose directly
ToothGrowth$dose = log2(ToothGrowth$dose)

# Classical analysis for comparison
lmToothGrowth <- lm(len ~ supp + dose + supp:dose, data=ToothGrowth)
summary(lmToothGrowth)

## ----GLMs----------------------------------------------------------------
full <- lmBF(len ~ supp + dose + supp:dose, data=ToothGrowth)
noInteraction <- lmBF(len ~ supp + dose, data=ToothGrowth)
onlyDose <- lmBF(len ~ dose, data=ToothGrowth)
onlySupp <- lmBF(len ~ supp, data=ToothGrowth)

allBFs <- c(full, noInteraction, onlyDose, onlySupp)
allBFs

## ----GLMs2---------------------------------------------------------------
full / noInteraction

## ----GLMposterior1-------------------------------------------------------
chainsFull <- posterior(full, iterations = 10000)

# summary of the "interesting" parameters
summary(chainsFull[,1:7])

## ----GLMposterior2,results='hide',echo=FALSE-----------------------------
chainsNoInt <- posterior(noInteraction, iterations = 10000)

## ----GLMplot,echo=FALSE,fig.width=10, fig.height=5-----------------------
ToothGrowth$dose <- ToothGrowth$dose - mean(ToothGrowth$dose)

cmeans <- colMeans(chainsFull)[1:6]
ints <- cmeans[1] + c(-1, 1) * cmeans[2]
slps <- cmeans[4] + c(-1, 1) * cmeans[5]


par(cex=1.8, mfrow=c(1,2))
plot(len ~ dose, data=ToothGrowth, pch=as.integer(ToothGrowth$supp)+20, bg = rgb(as.integer(ToothGrowth$supp)-1,2-as.integer(ToothGrowth$supp),0,.5),col=NULL,xaxt="n",ylab="Tooth length",xlab="Vitamin C dose (mg)")
abline(a=ints[1],b=slps[1],col=2)
abline(a=ints[2],b=slps[2],col=3)

axis(1,at=-1:1,lab=2^(-1:1))

dataVC <- ToothGrowth[ToothGrowth$supp=="VC",]
dataOJ <- ToothGrowth[ToothGrowth$supp=="OJ",]
lmVC <- lm(len ~ dose, data=dataVC)
lmOJ <- lm(len ~ dose, data=dataOJ)
abline(lmVC,col=2,lty=2)
abline(lmOJ,col=3,lty=2)

mtext("Interaction",3,.1,adj=1,cex=1.3)


# Do single slope

cmeans <- colMeans(chainsNoInt)[1:4]
ints <- cmeans[1] + c(-1, 1) * cmeans[2]
slps <- cmeans[4] 


plot(len ~ dose, data=ToothGrowth, pch=as.integer(ToothGrowth$supp)+20, bg = rgb(as.integer(ToothGrowth$supp)-1,2-as.integer(ToothGrowth$supp),0,.5),col=NULL,xaxt="n",ylab="Tooth length",xlab="Vitamin C dose (mg)")
abline(a=ints[1],b=slps,col=2)
abline(a=ints[2],b=slps,col=3)

axis(1,at=-1:1,lab=2^(-1:1))

mtext("No interaction",3,.1,adj=1,cex=1.3)

## ----eval=FALSE----------------------------------------------------------
#  chainsNoInt <- posterior(noInteraction, iterations = 10000)
#  
#  # summary of the "interesting" parameters
#  summary(chainsNoInt[,1:5])

## ----echo=FALSE----------------------------------------------------------
summary(chainsNoInt[,1:5])

## ------------------------------------------------------------------------
ToothGrowth$doseAsFactor <- factor(ToothGrowth$dose)
levels(ToothGrowth$doseAsFactor) <- c(.5,1,2)

aovBFs <- anovaBF(len ~ doseAsFactor + supp + doseAsFactor:supp, data = ToothGrowth)

## ------------------------------------------------------------------------
allBFs <- c(aovBFs, full, noInteraction, onlyDose)

## eliminate the supp-only model, since it performs so badly
allBFs <- allBFs[-1]

## Compare to best model
allBFs / max(allBFs)

## ----GLMplot2,echo=FALSE,fig.width=10, fig.height=5----------------------
plot(allBFs / max(allBFs))

## ------------------------------------------------------------------------
plot(Sepal.Width ~ Sepal.Length, data = iris)
abline(lm(Sepal.Width ~ Sepal.Length, data = iris), col = "red")

## ------------------------------------------------------------------------
cor.test(y = iris$Sepal.Length, x = iris$Sepal.Width)

## ------------------------------------------------------------------------
bf = correlationBF(y = iris$Sepal.Length, x = iris$Sepal.Width)
bf

## ------------------------------------------------------------------------
samples = posterior(bf, iterations = 10000)

## ------------------------------------------------------------------------
summary(samples)

## ------------------------------------------------------------------------
plot(samples[,"rho"])

## ----propprior,echo=FALSE,fig.width=10, fig.height=5---------------------
p0 = .5
rnames = c("medium","wide","ultrawide")
r = sapply(rnames,function(rname) BayesFactor:::rpriorValues("proptest",,rname))
leg_names = paste(rnames," (r=",round(r,3), ")", sep="")

omega = seq(-5,5,len=100)
pp = dlogis(omega,qlogis(p0),r[1])

plot(omega,pp, col="black", typ = 'l', lty=1, lwd=2, ylab="Prior density", xlab=expression(paste("True log odds ", omega)), yaxt='n')

pp = dlogis(omega,qlogis(p0),r[2])
lines(omega, pp, col = "red",lty=1, lwd=2)

pp = dlogis(omega,qlogis(p0),r[3])
lines(omega, pp, col = "blue",lty=1,lwd=2)

axis(3,at = -2:2 * 2, labels=round(plogis(-2:2*2),2))
mtext(expression(paste("True probability ", pi)),3,2,adj=.5)

legend(-5,.5,legend = leg_names, col=c("black","red","blue"), lwd=2,lty=1)

## ------------------------------------------------------------------------
bf = proportionBF( 682, 682 + 243, p = 3/4)
1 / bf

## ------------------------------------------------------------------------
binom.test(682, 682 + 243, p = 3/4)

## ----proppost,fig.width=10, fig.height=5---------------------------------
chains = posterior(bf, iterations = 10000)
plot(chains[,"p"], main = "Posterior of true probability\nof 'giant' progeny")

## ----results='asis', echo=FALSE------------------------------------------
data(raceDolls)
kable(raceDolls)

## ------------------------------------------------------------------------
bf = contingencyTableBF(raceDolls, sampleType = "indepMulti", fixedMargin = "cols")
bf

## ------------------------------------------------------------------------
chisq.test(raceDolls)

## ------------------------------------------------------------------------
chains = posterior(bf, iterations = 10000)

## ------------------------------------------------------------------------
sameRaceGivenWhite = chains[,"pi[1,1]"] / chains[,"pi[*,1]"]
sameRaceGivenBlack = chains[,"pi[1,2]"] / chains[,"pi[*,2]"]


## ----ctablechains,fig.width=10, fig.height=5-----------------------------
plot(mcmc(sameRaceGivenWhite - sameRaceGivenBlack), main = "Increase in probability of child picking\nsame race doll (white - black)")

## ------------------------------------------------------------------------
data(puzzles)

puzzleGenBF <- generalTestBF(RT ~ shape + color + shape:color + ID, data=puzzles, whichRandom="ID")

puzzleGenBF

## ------------------------------------------------------------------------
puzzleGenBF <- generalTestBF(RT ~ shape + color + shape:color + ID, data=puzzles, whichRandom="ID", neverExclude="ID")

puzzleGenBF

## ------------------------------------------------------------------------
puzzleGenBF <- generalTestBF(RT ~ shape + color + shape:color + shape:ID + ID, data=puzzles, whichRandom="ID", neverExclude="ID")

puzzleGenBF

## ------------------------------------------------------------------------
puzzleGenBF <- generalTestBF(RT ~ shape + color + shape:color + shape:ID + ID, data=puzzles, whichRandom="ID", neverExclude="^ID$")

puzzleGenBF

## ------------------------------------------------------------------------
puzzleCullBF <- generalTestBF(RT ~ shape + color + shape:color + ID, data=puzzles, whichRandom="ID", noSample=TRUE,whichModels='all')

puzzleCullBF

## ------------------------------------------------------------------------
missing = puzzleCullBF[ is.na(puzzleCullBF) ]
done = puzzleCullBF[ !is.na(puzzleCullBF) ]

missing

## ------------------------------------------------------------------------
# get the names of the numerator models
missingModels = names(missing)$numerator

# search them to see if they contain "shape" or "color" - 
# results are logical vectors
containsShape = grepl("shape",missingModels)
containsColor = grepl("color",missingModels)

# anything that does not contain "shape" and "color"
containsOnlyOne = !(containsShape & containsColor)

# restrict missing to only those of interest
missingOfInterest = missing[containsOnlyOne]

missingOfInterest

## ------------------------------------------------------------------------
# recompute the Bayes factors for the missing models of interest
sampledBayesFactors = recompute(missingOfInterest)

sampledBayesFactors

# Add them together with our other Bayes factors, already computed:
completeBayesFactors = c(done, sampledBayesFactors)

completeBayesFactors

## ------------------------------------------------------------------------
data(puzzles)

# Get MCMC chains corresponding to "full" model
# We prevent sampling so we can see the parameter names
# iterations argument is necessary, but not used
fullModel = lmBF(RT ~ shape + color + shape:color + ID, data = puzzles, noSample=TRUE, posterior = TRUE, iterations=3)

fullModel

## ------------------------------------------------------------------------
fullModelFiltered = lmBF(RT ~ shape + color + shape:color + ID, data = puzzles, noSample=TRUE, posterior = TRUE, iterations=3,columnFilter="ID")

fullModelFiltered

## ------------------------------------------------------------------------
# Sample 10000 iterations, eliminating ID columns
chains = lmBF(RT ~ shape + color + shape:color + ID, data = puzzles, posterior = TRUE, iterations=10000,columnFilter="ID")

## ----acfplot,fig.width=10,fig.height=5,echo=FALSE------------------------
par(mfrow=c(1,2))
plot(as.vector(chains[1:1000,"shape-round"]),type="l",xlab="Iterations",ylab="parameter shape-round")
acf(chains[,"shape-round"])

## ------------------------------------------------------------------------
chainsThinned = recompute(chains, iterations=20000, thin=2)

# check size of MCMC chain
dim(chainsThinned)

## ----acfplot2,fig.width=10,fig.height=5,echo=FALSE-----------------------
par(mfrow=c(1,2))
plot(as.vector(chainsThinned[1:1000,"shape-round"]),type="l",xlab="Iterations",ylab="parameter shape-round")
acf(chainsThinned[,"shape-round"])

## ----tidy=FALSE----------------------------------------------------------
newprior.bf = anovaBF(RT ~ shape + color + shape:color + ID, data = puzzles,
                           whichRandom = "ID",rscaleEffects = c( color = 1 ))

newprior.bf

