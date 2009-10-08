dG <- function(x,data,prior=1,logg=1)
{
	.Call("RdG",x,data,as.integer(prior),as.integer(logg),package="sampG")
}

