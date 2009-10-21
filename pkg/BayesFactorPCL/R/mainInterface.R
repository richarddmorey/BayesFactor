
BFAnovaGUI <- function(data=NULL,filename=NULL,setup=NULL,skipInstr=FALSE,name="Analysis",interface=TRUE,doOutput=TRUE, storePred=FALSE){

  noMoreTesting=0
  BFAnovaInit()
  passedPackage = class(setup)=="PCLBFAnova"
  pack = new("PCLBFAnova")
 
  if(!(skipInstr | !is.null(setup) | !interface)){ 
    a=BFAnovaInitInstructions()
  }else{
    a=0
  }
  if(is.null(a)) {cat("Cancelled by user.\n");return(1)}
  
  # Load in data
  if(is.null(data) & is.null(setup)){
    if(!is.null(filename)){
      data=read.csv(filename)
    }else{
      filename <- tclvalue(tkgetOpenFile(
        filetypes = "{{CSV files} {.csv}} {{All files} *}"))
    	if (filename == "") stop("No file selected.");
    	data <- read.csv(filename)
    }
  }
     
  if(passedPackage){
    pack@settings = setup@settings
    rm(setup)
    pack@settings@time = date()
    pack@settings@filename = "Settings loaded from PCLBFAnova object."
      
  }else{
    pack@settings=new("PCLBFAnovaSettings")
    
    pack@settings@data         = data
    pack@settings@filename     = as.character(filename)
    pack@settings@analysisname = as.character(name)
    pack@settings@time         = date()
    
	 
    pack@settings@SelCols  = BFAnovaSelectColumns(pack@settings@data)
      if(identical(pack@settings@SelCols,list(0)))  {cat("Cancelled by user.\n");return(1)}
	
    
    pack@settings@newDat      = createDat(pack@settings@data,pack@settings@SelCols$dependent,pack@settings@SelCols$selectedcols)
    pack@settings@mods        = listModels(pack@settings@newDat,pack@settings@SelCols$selectedcols)
    pack@settings@effs        = niceListEffects(pack@settings@mods)
     
    pack@settings@SelEffs     = BFAnovaSelectEffects(pack@settings@effs)
      if(identical(pack@settings@SelEffs,list(0))) {cat("Cancelled by user.\n");return(1)}
    
  
    pack@settings@intMods     = pack@settings@SelEffs[[3]][pack@settings@SelEffs[[3]][,5]==TRUE,1:2]
    pack@settings@SelEffs2    = pack@settings@SelEffs[[3]][pack@settings@SelEffs[[3]][,5]==TRUE,]
    pack@settings@newDat2  = createModCols(pack@settings@newDat,pack@settings@mods,pack@settings@intMods,pack@settings@SelCols$selectedcols)

   pack@settings@designMatrix = createDesignMatrix(pack@settings@newDat2) 
   pack@settings@namedDat2   = createMeaningfulCols(pack@settings@newDat,pack@settings@mods,pack@settings@intMods,pack@settings@SelCols$selectedcols)
    
    pack@settings@Lvls        =  apply(pack@settings@newDat2Cat[,-(1:5)],2,function(v) length(unique(v)))
                             
    pack@settings@PriorSetup = BFAnovaPriorSetup(useA=pack@settings@useA,minWishDF=minWishDF)	
	if(identical(pack@settings@PriorSetup,list(0))) {cat("Cancelled by user.\n");return(1)}

    
    }      
      if(interface){
        if(!identical(pack@settings@MCMCSetup,list())){
          pack@settings@MCMCSetup = BFAnovaMCMCSetup(nIters=pack@settings@MCMCSetup$nIter,
            burnin=pack@settings@MCMCSetup$burnin,
            progress=pack@settings@MCMCSetup$progress)
        }else{
          pack@settings@MCMCSetup = BFAnovaMCMCSetup(nIters=1000,burnin=200,progress=10)
        }
      }
      if(identical(pack@settings@MCMCSetup,list(0))) {cat("Cancelled by user.\n");return(1)}
      
      
      if(interface & doOutput)
        {
            OutputSetup = BFAnovaOutputSetup()
            if(identical(OutputSetup,list(0))) {cat("Cancelled by user.\n");return(1)}
			doOutput=as.integer(OutputSetup$doOutput)
			storePred = as.integer(OutputSetup$Pred)
		}
		
     
      ## Create progress bar
      progress=as.integer(pack@settings@MCMCSetup$progress)
      if(progress){ pb = txtProgressBar(min = 0, max = as.integer(pack@settings@MCMCSetup$nIter), style = 3) }else{ pb=NULL }
      pbFun = function(samps){ if(progress) setTxtProgressBar(pb, samps)}
      

      pack@output = new("PCLBFAnovaOutput")
      
      cat("Starting MCMC...\n")
      

        output = .Call("WM2_GibbsSamplerNoCov", ..., PACKAGE="")
      
      pack@output@burnin = as.numeric(pack@settings@MCMCSetup$burnin)
      pack@output@par = #

      chainnames=paste(pack@output@par[,4],pack@output@par[,2],pack@output@par[,3],sep=" ")
      chainnames=paste(chainnames,pack@output@par[,1],sep=" on ")
      pack@output@Effchains=mcmc(t(output[[1]]))
      dimnames(pack@output@Effchains)[[2]] = chainnames
      


  if(doOutput){
    cat("Writing output files...\n")
    flush.console()
    writeOutputFiles(pack,directory=getwd(),outputSettings=OutputSetup,storePred)
  }
  return(pack)
}

