
#'Opens the HTML manual for the BayesFactor package
#'
#'This function opens the HTML manual for the BayesFactor package in whatever
#'browser is configured.
#'
#'This function opens the HTML manual for the BayesFactor package in whatever
#'browser is configured.
#'@return \code{BFManual} returns \code{NULL} invisibly.
#'@author Richard D. Morey (\email{richarddmorey@@gmail.com})
#'@keywords misc
#'@export
BFManual <- function(){
  myPort <- tools:::httpdPort
  
  # myPort will be 0 if dynamic help has not been started
  if(!myPort) myPort <- tools::startDynamicHelp()
  
  browseURL(paste0('http://127.0.0.1:', myPort,'/library/BayesFactor/doc/index.html'))
  invisible(NULL)
}

