\name{nWayAOV2}
\alias{nWayAOV2}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Backend function to compute the Bayes factors in ANOVA designs
}
\description{
This function computes Bayes factors. It is a utility function not designed for end-users.
}
\usage{
nWayAOV2(modNum, env, rscaleFixed, rscaleRandom, progressCallback=NULL, ...)
    
}
\arguments{
%- maybe also 'usage' for other objects documented here.
\item{modNum}{
An integer model number specifying which model to compute the Bayes factor for
}
\item{env}{
An environment containing the data and design matrices for the model fit
}
\item{progressCallback}{

}
\item{rscaleFixed}{The prior scale of effect size for factors that are treated as a fixed.
}
\item{rscaleRandom}{
The prior scale of effect size for factors that are treated as random.}

\item{...}{
Additional arguments to pass to \code{\link{nWayAOV.MC}}
}
}

\details{
This model is a utility function for the BayesFactor package. It is not designed to called by end-users.
}


\value{
This function returns a Bayes factor against the null model in which all fixed effects are null.
}

\author{
Richard D. Morey
}



%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{allNways}}, \code{\link{nWayAOV.MC}}, \code{\link{aov}}, \code{\link{ttest.Quad}}, \code{\link{oneWayAOV.Quad}},\code{\link{oneWayAOV.Gibbs}}, \code{\link[doMC:doMC-package]{doMC}}, \code{\link{options}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
