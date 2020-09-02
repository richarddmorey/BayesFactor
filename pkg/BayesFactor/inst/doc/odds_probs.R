## ----echo=FALSE,message=FALSE,results='hide'-----------------------------
options(markdown.HTML.stylesheet = 'extra/manual.css')
library(knitr)
options(digits=3)
require(graphics)
set.seed(2)

## ----message=FALSE,results='hide',echo=FALSE-----------------------------
library(BayesFactor)
options(BFprogress = FALSE)
bfversion = BFInfo()
session = sessionInfo()[[1]]
rversion = paste(session$version.string," on ",session$platform,sep="")

## ------------------------------------------------------------------------
data(puzzles)
bf = anovaBF(RT ~ shape*color + ID, whichRandom = "ID", data = puzzles)
bf

## ------------------------------------------------------------------------
prior.odds = newPriorOdds(bf, type = "equal")
prior.odds

## ------------------------------------------------------------------------
priorOdds(prior.odds) <- c(4,3,2,1)
prior.odds

## ------------------------------------------------------------------------
post.odds = prior.odds * bf
post.odds

## ------------------------------------------------------------------------
post.prob = as.BFprobability(post.odds)
post.prob

## ------------------------------------------------------------------------
post.prob / .5

## ------------------------------------------------------------------------
post.prob[1:3]

## ------------------------------------------------------------------------
post.prob[1:3] / 1

