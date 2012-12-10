\name{puzzles}
\alias{puzzles}
\docType{data}
\title{Puzzle completion times from Hays (1994)}
\description{
	Puzzle completion time example data from Hays (1994).
}
\usage{data(puzzles)}
\format{
  A data frame with 48 observations on 3 variables.
  \describe{
  	\item{\code{RT}}{Puzzle completion time, in minutes}
    \item{\code{ID}}{the subject identifier}
    \item{\code{shape}}{shape of the puzzle (round or square)}
    \item{\code{color}}{color content of the puzzle (monochromatic or color)}
  }
}
\details{
	Hays (1994; section 13.21, table 13.21.2, p. 570) describes a experiment wherein 12 participants complete four puzzles each. Puzzles could be either square or round, and either monochromatic or in color. Each participant completed every combination of the two factors.
}
\source{
Hays, W. L. (1994), Statistics (5th edition), Harcourt Brace, Fort Worth, Texas
}
\examples{
data(puzzles)

## classical ANOVA
## Both color and shape are significant, interaction is not
classical <- aov(RT ~ shape*color + Error(ID/(shape*color)), data=puzzles)
summary(classical)

## Bayes Factor
## Best model is main effects model, no interaction
allNways(y = puzzles$RT, dataFixed = puzzles[,3:4], dataRandom = puzzles$ID, progress = FALSE)

}
\keyword{datasets}
