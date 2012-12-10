\name{aovGUI}
\alias{aovGUI}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to start GUI for computing the Bayes factors in an ANOVA design
}
\description{
This function starts a graphical user interface that allows model comparisons with data from ANOVA designs. 
}
\usage{

aovGUI(y, dataFixed=NULL, dataRandom=NULL, iterations = 1000, 
    		rscaleFixed=.5, rscaleRandom=1)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	A vector of observations	
}
  \item{dataFixed}{
A data frame containing the fixed factors, one column per factor. Columns will be converted to factors if they are not factors already.
}
  \item{dataRandom}{
A data frame containing the random factors, one column per factor. Columns will be converted to factors if they are not factors already.
}
  \item{iterations}{
Number of Monte Carlo samples to perform (see details).  
}
\item{rscaleFixed}{The prior scale of effect size for factors that are treated as a fixed.
}
\item{rscaleRandom}{
The prior scale of effect size for factors that are treated as random.}
}

\details{
	Thus function uses the Rook webserver in R to build a user interface. It will open the interface in a new browser window. Instructions for the user interface can be found by clicking 'Help' once the interface starts.
	
	See the help for \code{\link{allNways}} and \code{\link{nWayAOV.MC}} for more details.
}


\value{
\code{aovGUI} returns nothing.
}

\author{
Richard D. Morey
}



%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{
This interface is still in a very early stage of development.
}
\seealso{
\code{\link{allNways}},\code{\link{nWayAOV.MC}}, \code{\link{aov}}
}
\examples{
data(puzzles)

\dontrun{
aovGUI(y = puzzles$RT, dataFixed = puzzles[,3:4], dataRandom = puzzles$ID)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
