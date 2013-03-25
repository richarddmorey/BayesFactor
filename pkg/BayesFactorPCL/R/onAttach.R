.onAttach<- function(libname, pkgname){
  packageStartupMessage("************\nWelcome to ",pkgname," ",BFInfo(FALSE),". If you have",
                        " questions, please contact Richard Morey (richarddmorey@gmail.com).\n\n",
                        "Type BFManual() to open the manual.\n************", 
                        appendLF = TRUE)
  setOptions()
}

setOptions <- function(){
  
  if(is.null(options()$BFMaxModels)) options(BFMaxModels = 50000)
  if(is.null(options()$BFpretestIterations)) options(BFpretestIterations = 100)
  if(is.null(options()$BFapproxOptimizer)) options(BFapproxOptimizer = "optim")
  if(is.null(options()$BFapproxLimits)) options(BFapproxLimits = c(-15,15))
  if(is.null(options()$BFprogress)) options(BFprogress = TRUE)
  
}

