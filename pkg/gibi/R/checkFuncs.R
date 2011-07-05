.checkNonNegative(v,...)
	{
		if(.checkNumConvert(v)){
			return(all(as.numeric(v) >= 0))
		}else{
			return(FALSE)
		}
	}

.checkPositive(v,...)
	{
		if(.checkNumConvert(v)){
			return(all(as.numeric(v) > 0))
		}else{
			return(FALSE)
		}
	}

.checkConvert <- function(v, type, checkFunc)
	{
		saveOpt = options()$warn
		options(warn=2)
		num = try(checkFunc(v),silent=TRUE)
		options(warn=saveOpt)
		return(class(num) == type)
	}

.checkCharConvert <- function(v,...)
	{
		.checkConvert(v, "character", as.character)	
	}

.checkIntConvert <- function(v,...)
	{
		.checkConvert(v, "integer", as.integer)	
	}

.checkNumConvert <- function(v,...)
	{
		.checkConvert(v, "numeric", as.numeric)
	}

.checkFacConvert <- function(v,...)
	{
		.checkConvert(v, "factor", as.factor)
	}

.checkNumLevels <- function(v,minLevs=NULL,maxLevs=NULL,reqLevs=NULL)
	{
		if(.checkFacConvert(v)){
			retVal=TRUE
			nlevs = nlevels(as.factor(v))
			if(!is.null(maxLevs)){
				retVal = retVal & (nLevs <= maxLevs)
			}
			if(!is.null(minLevs)){
				retVal = retVal & (nLevs >= minLevs)
			}
			if(!is.null(reqLevs)){
				retVal = retVal & (nLevs %in% reqLevs)
			}
			if(is.null(minLevs) & is.null(maxLevs) & is.null(reqLevs)){
				warning("In .checkNumLevels(), no criteria set; returning TRUE.")
			}
			return(retVal)
		}else{
			return(FALSE)
		}
	}
