#'Prints the version information for the BayesFactor package
#'
#'Prints the version, revision, and date information for the BayesFactor package
#'
#'This function prints the version and revision information for the BayesFactor
#'package.
#'
#'@param print if \code{TRUE}, print version information to the console
#'@return \code{BFInfo} returns a character string containing the version and
#'  revision number of the package..
#'@author Richard D. Morey (\email{richarddmorey@@gmail.com})
#'@keywords misc
#'@export
BFInfo <- function(print=TRUE)
{
  fn <- system.file("SVN_VERSION", package="BayesFactor") 
  if (file.exists(fn)) { 
    svn_version <- scan(system.file("SVN_VERSION", package="BayesFactor"), 
                        what=character(1), sep="\n", quiet=TRUE) 
  } else { 
    svn_version <- "(unknown)" 
  }
  if(print){
    cat("Package BayesFactor\n")
	  cat(packageDescription("BayesFactor")$Version,"\n")
    cat("SVN revision:", svn_version,"\n")
  }
  retStr = paste(packageDescription("BayesFactor")$Version, "//", svn_version)
  invisible(retStr)
} 
