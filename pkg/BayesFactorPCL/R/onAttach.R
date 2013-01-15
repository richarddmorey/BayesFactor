.onAttach<- function(libname, pkgname){
  packageStartupMessage("************\nWelcome to ",pkgname," ",BFInfo(FALSE),". Much has changed since the ",
                        "last version, including new functionality and function name changes. ",
                        "Please see the help files for the new function names and arguments. If you have",
                        " questions, please contact Richard Morey (richarddmorey@gmail.com).\n\n",
                        "Type BFManual() to open the manual.\n************", 
                        appendLF = TRUE)
}