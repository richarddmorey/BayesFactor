BFInfo <- function()
{
  fn <- system.file("SVN_VERSION", package="BayesFactor") 
  if (file.exists(fn)) { 
    svn_version <- scan(system.file("SVN_VERSION", package="BayesFactor"), 
                        what=character(1), sep="\n", quiet=TRUE) 
  } else { 
    svn_version <- "(unknown)" 
  } 
  cat("Package BayesFactor\n")
	cat(packageDescription("BayesFactor")$Version,"\n")
  cat("SVN revision:", svn_version,"\n")
  invisible()
} 