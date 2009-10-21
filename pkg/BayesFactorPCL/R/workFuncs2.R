
int.matrix<-function(m){
  dims=dim(m)
  m2=as.integer(as.matrix(m))
  if(is.null(dims)){
    dim(m2)=c(length(m),1)
  }else{  
    dim(m2)=dims
  }
  return(m2)
}



createDat <- function(data,depColumn,interestColumns)
{

	names = colnames(data)
	depColumnNum = match(depColumn,names)
	interestColumnsNum = match(interestColumns[,1],names)

	newData = data[,c(depColumnNum,interestColumnsNum)]
	
	return(newData)
}


factMax = function(v){
	if(is.integer(v)){
		max(v)
	}else{
		1
	}
}

listModels <- function(newDat,interestColumns){
	nFactors=dim(newDat)[2]-1
	names=colnames(newDat)[-1]
	levels=1:nFactors
	for(i in 1:nFactors){
	  levels[i]=length(unique(newDat[,5+i])) 
	}
	modsInt=list()

	for(i in 1:nFactors){
		modsInt[[i]]=combinations(nFactors,i)
	}

	list(modsInt,names,levels)
}

niceListEffects = function(mods){
	ret=NULL
	for(i in 1:length(mods[[1]])){
		for(j in 1:dim(mods[[1]][[i]])[1]){
			cols=mods[[1]][[i]][j,]
			if(length(cols)==1){
				ret=rbind(ret,c(i,j,mods[[3]][mods[[1]][[i]][j]],mods[[2]][mods[[1]][[i]][j]]))
			}else{
				ret=rbind(ret,c(i,j,prod(mods[[3]][mods[[1]][[i]][j,]]),paste(mods[[2]][mods[[1]][[i]][j,]],collapse=" by ")))
			}

		}
	}
	ret=data.frame(ret)
	ret[,1]=as.numeric(as.character(ret[,1]))
	ret[,2]=as.numeric(as.character(ret[,2]))
	ret[,3]=as.numeric(as.character(ret[,3]))
	ret=data.frame(ret,ret[,1]*0)
	colnames(ret)=c("way","effnum","nlevels","name","S")
	ret
}


createModCols <- function(newDat,allMods,intMods,SelCols){

	if(is.null(dim(intMods))){
	  nNewCols = length(intMods)
	}else{
	  nNewCols = dim(intMods)[1]
	}
	#nNewCols = dim(intMods)[1]
	newDat2  = newDat[,1]

	names=c("Dependent")

	for(i in 1:nNewCols){
      		myCols=allMods[[1]][[intMods[i,1]]][intMods[i,2],]
      		nLvls=allMods[[3]][myCols]
      		if(length(myCols)>1){
			newDat2=data.frame(newDat2,as.integer(as.factor(apply(cbind(newDat[,(myCols+1)]),1,paste,collapse="1"))))
	      		names=c(names,paste(colnames(newDat)[(myCols+1)[cc=="CATEG"]],collapse=".x."))				      			
	      		}
		}else{

			newDat2=data.frame(newDat2,as.integer(as.factor(newDat[,myCols+1])))
      			names=c(names,colnames(newDat)[myCols+1])

      		}
	}
colnames(newDat2)=names
return(newDat2)
}


createMeaningfulCols <- function(newDat,allMods,intMods,SelCols)
{

	if(is.null(dim(intMods))){
	  nNewCols = length(intMods)
	}else{
	  nNewCols = dim(intMods)[1]
	}

	newDat2  = newDat[,1]
	namedCols=data.frame(newDat[,-1])


	names=c("Dependent")

	for(i in 1:nNewCols){
      		myCols=allMods[[1]][[intMods[i,1]]][intMods[i,2],]
      		nLvls=allMods[[3]][myCols]
      		cc=CatOrCont(colnames(newDat)[myCols+5],SelCols)
      		if(length(myCols)>1){
			newDat2=data.frame(newDat2,apply(cbind(namedCols[,myCols]),1,paste,collapse=".x."))
			names=c(names,paste(colnames(newDat)[myCols+1],collapse=".x.")) 
		}else{
			newDat2=data.frame(newDat2,namedCols[,myCols])
      			names=c(names,colnames(newDat)[myCols+1]) 
      		}
	}
colnames(newDat2)=names
return(newDat2)
}


createDesignMatrix(newDat2){
	nCols = dim(newDat2)[2]-1
	nRows = dim(newDat2)[1]
	X = matrix(1,nrow=nRows,ncol=1)
	for(i in 1:nCols){
	      nLevs = length(unique(newDat2[,i+1]))
	      thisCol = newDat2[,i+1]
	      X0=matrix(0,nrow=nRows,ncols=nLevs)
	      X0[nRows*(thisCol-1)+0:(nRows-1)+1]=1
	      X=cbind(X,X0)
	}
return(X)
}

getLevels = function(data,SelCols)
{
	colnums=match(SelCols[,1],colnames(data))
	myLevels=list()
	for(i in 1:length(colnums))
	{
		if(SelCols[i,2]=="CATEG"){
			myLevels[[i]]=levels(as.factor(data[,colnums[i]]))
		}else{
			myLevels[[i]]=""
		}
	}
myLevels
}

getNFactorLevels <-function(allMods,intMods)
{
	nNewCols = dim(intMods)[1]

	FactorLevels=array(1:nNewCols*NA)
	names=NULL

	for(i in 1:nNewCols){
      		myMod=allMods[[1]][[intMods[i,1]]][intMods[i,2],]
      		if(length(myMod)>1){
			FactorLevels[i]=prod(allMods[[3]][myMod])
      			names=c(names,paste(allMods[[2]][myMod],collapse=".x."))
	      	}else{
			FactorLevels[i]=allMods[[3]][myMod]
      			names=c(names,allMods[[2]][myMod])
      		}
	}

	names(FactorLevels)=names
	FactorLevels
}


