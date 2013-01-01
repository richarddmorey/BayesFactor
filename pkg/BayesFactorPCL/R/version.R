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