rG <-
function(start=0,N,data,prior=1,method=1,logg=1,control=1,progress=0)
{
	## Create progress bar
    if(progress){ pb = txtProgressBar(min = 0, max = as.integer(N), style = 3) }else{ pb=NULL }
    pbFun = function(samps){ if(progress) setTxtProgressBar(pb, samps)}
      
	ret = .Call("RsampG",start,as.integer(N),data,as.integer(prior),as.integer(method),
	        as.integer(logg), control, progress, pbFun, new.env(), package="sampG")
	
	cat("\n")
	return(ret)
}

