.onAttach<- function(libname, pkgname){
  packageStartupMessage("************\nWelcome to ",pkgname," ",BFInfo(FALSE),". Much has changed since the ",
                        "last version, including new functionality and function name changes. ",
                        "Please see the help files for the new function names and arguments.\n************", 
                        appendLF = TRUE)
}