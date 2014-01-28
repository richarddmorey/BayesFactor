.onAttach<- function(libname, pkgname){
  packageStartupMessage("************\nWelcome to ",pkgname," ",BFInfo(FALSE),". If you have",
                        " questions, please contact Richard Morey (richarddmorey@gmail.com).\n\n",
                        "Type BFManual() to open the manual.\n************", 
                        appendLF = TRUE)
  setOptions()
}

#'options() for package BayesFactor
#'
#'Options that can be set for the BayesFactor package
#'
#'The BayesFactor package has numerous options that can be set to globally 
#'change the behavior of the functions in the package. These options can be
#'changed using \code{\link[base]{options}}().
#'
#'\describe{
#'\item{\code{BFMaxModels}}{Integer; maximum number of models to analyze in \code{\link{anovaBF}} or \code{\link{regressionBF}}} 
#'\item{\code{BFprogress}}{If \code{TRUE}, progress bars are on by default; if \code{FALSE}, they are disabled by default.}
#'\item{\code{BFpretestIterations}}{Integer; if sampling is needed to compute the Bayes factor, the package attempts to 
#'choose the most efficient sampler. This option controls the number of initial test iterations.}
#'\item{\code{BFapproxOptimizer}}{\code{"nlm"} or \code{"optim"}; changes the optimization function used for the importance sampler. If one fails, try the other.}
#'\item{\code{BFapproxLimits}}{Vector of length two containing the lower and upper limits on 
#'on \code{log(g)} before the the posterior returns \code{-Inf}. This only affects the initial optimization step for the importance sampler.}
#'\item{\code{BFfactorsMax}}{Maximum number of factors to try to do enumeration with in generalTestBF.}
#'}
#'
#'@name options-BayesFactor
#'@seealso \code{\link[base]{options}}

NULL



setOptions <- function(){
  
  if(is.null(options()$BFMaxModels)) options(BFMaxModels = 50000)
  if(is.null(options()$BFpretestIterations)) options(BFpretestIterations = 100)
  if(is.null(options()$BFapproxOptimizer)) options(BFapproxOptimizer = "optim")
  if(is.null(options()$BFapproxLimits)) options(BFapproxLimits = c(-15,15))
  if(is.null(options()$BFprogress)) options(BFprogress = TRUE)
  if(is.null(options()$BFfactorsMax)) options(BFfactorsMax = 5)
  
}

