\name{jnt.post}
\Rdversion{1.1}
\alias{jnt.post}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
jnt.post(eta2, delta, t, N, df.prior = 1, r.scale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta2}{
%%     ~~Describe \code{eta2} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{df.prior}{
%%     ~~Describe \code{df.prior} here~~
}
  \item{r.scale}{
%%     ~~Describe \code{r.scale} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (eta2, delta, t, N, df.prior = 1, r.scale = 1) 
exp(-(N/2 + 1) * log(eta2) - (N - 1)/(2 * eta2) - (t - delta * 
    sqrt(eta2 * N))^2/(2 * eta2) - (df.prior + 1)/2 * log(1 + 
    (delta/r.scale)^2/df.prior))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
