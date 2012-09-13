\name{logMeanExpLogs}
\alias{logMeanExpLogs}
\alias{logCumMeanExpLogs}

\title{
Functions to compute the logarithm of the mean (and cumulative means) of vectors of logarithms 
}
\description{
Given a vector of numeric values of real values represented in log form, \code{logMeanExpLogs} computes the logarithm of the mean of the (exponentiated) values. \code{logCumMeanExpLogs} computes the logarithm of the cumulative mean.

}
\usage{
logMeanExpLogs(v)

logCumMeanExpLogs(v)
}
\arguments{
  \item{v}{
A vector of (log) values
}
}
\details{
Given a vector of values of log values \var{v}, one could compute \code{log(mean(exp(v)))} in R. However, exponentiating and summing will cause a loss of precision, and possibly an overflow. These functions use the identity
\deqn{\log(e^a + e^b) = a + \log(1+e^{b-a})}{%
     log(e^a + e^b) = a + log[ 1 + e^(b-a) ]}
and the method of computing \eqn{\log(1+e^x)}{log(1+e^x} that avoids overflow (see the references). The code is written in C for very fast computations.
}
\value{
\code{logMeanExpLogs} returns a single value; \code{logCumMeanExpLogs} returns a vector of values of the same length as \var{v}.
}
\references{
For details of the approximation of \eqn{\log(1+e^x)}{log(1+e^x)} used to prevent loss of precision, see \url{http://www.codeproject.com/Articles/25294/Avoiding-Overflow-Underflow-and-Loss-of-Precision}.
}
\author{
Richard D. Morey
}
\examples{
# Sample 100 values
y = log(rexp(100,1))

# These will give the same value, 
# since e^y is "small"
logMeanExpLogs(y)
log(mean(exp(y)))

# We can make e^x overflow by multiplying
# e^y by e^1000
largeVals = y + 1000

# This will return 1000 + log(mean(exp(y)))
logMeanExpLogs(largeVals)

# This will overflow
log(mean(exp(largeVals)))
}
\keyword{arith}
\keyword{misc}