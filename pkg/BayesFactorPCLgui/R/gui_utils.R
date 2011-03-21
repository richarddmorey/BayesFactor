bfPCL_extractEnv <- function()
return(bfEnv)


.notebookPages <- function(name)
{
COLUMNS <- c(intro=0,data=1,analysis=2,settings=3,log=4)
as.integer(COLUMNS[name])
}


.fileSafeString <- function(string)
{
	gsub("[^\\w-\\s]","_",string,perl=TRUE)
}

theWidget <- function(name) {
	return(StateEnv$GUI$getObject(name))
}


freezeGUI <- function(echo.to.log=T) {
	StateEnv$win$setSensitive(F)
	StateEnv$win$getWindow()$setCursor(gdkCursorNew("watch"))
	#StateEnv$echo.to.log <- echo.to.log
	#setStatusBar("")
}


thawGUI <- function() {
	StateEnv$win$setSensitive(T)
	StateEnv$win$getWindow()$setCursor(NULL)
	#StateEnv$echo.to.log <- T # default
}

getpackagefile <- function(filename) {
	## Try firstly to load from the installed package
	## Otherwise, look locally.
	myPath <- system.file("etc", filename, package = "BayesFactorPCLgui")
	if (identical(myPath, "")) 
		myPath <- file.path("BayesFactorPCLgui", "BayesFactorPCLgui", "inst", 
			"etc", filename)
	if (!file.exists(myPath)) stop("could not find file ", filename)
	myPath
}

.setStatusBarText <- function(text1=NULL,context="General")
{
	status = theWidget("statusbar1")
	context.id=gtkStatusbarGetContextId(status,context)
	if(!is.null(text1)){
		gtkStatusbarPop(status, context.id)
		gtkStatusbarPush(status, context.id, text1)
	}
}

.present_main_window_after_destroy<-function(widget)
{
	StateEnv$win$present()
}

.saneNum <- function(x,prec=0)
{
	if(is.null(x)) return(NULL)
	as.character(round(x,prec))
}

.getDataFrameList<-function(envir=globalenv())
{
	n =  sapply(ls(env=envir),function(x){
			eval(parse(text=paste("is.data.frame(",x,")",sep="")),env=globalenv())
				}
	)
	
	if(length(n)==0){
		return(data.frame())
	}
	n=n[n]
	data.frame(Name=names(n))
}


.selectDataFrame <- function(envir = .GlobalEnv) 
{
	listOfObjects <- .getDataFrameList(envir=envir)
	if(dim(listOfObjects)[1]<1){
		gmessage(paste("There are no data frames in the R global environment."), title="Data error",
		icon="error",toolkit=guiToolkit("RGtk2"))
		return(0)

	}
	toplevel = gwindow("Select R data.frame...")
	gtable(listOfObjects,
		container = toplevel,
		action=NULL,
		handler = function(h,...) {
			.setDataForColumnSelection(eval(parse(text=as.character(svalue(h$obj)))))
			dispose(toplevel)
		},toplevel=toplevel)
}

.getDataFrame <- function(name)
{
	if(length(name)>1) return(NULL)
	z = eval(parse(text=name),env=globalenv())
	if(is.data.frame(z)){
		return(z)
	}else{
		return(NULL)
	}
}

.integerMatrix<-function(m){
  dims=dim(m)
  m2=as.integer(as.matrix(m))
  if(is.null(dims)){
    dim(m2)=c(length(m),1)
  }else{  
    dim(m2)=dims
  }
  return(m2)
}

.treeModelGetNthCol<-function(treemodel,n=0)
{
	totalIters = gtkTreeModelIterNChildren(treemodel, NULL)
	if(totalIters==0) return(NULL)


	retval = 1:totalIters * NA
	
	for(i in 1:totalIters)
	{
		iter = gtkTreeModelGetIterFromString(treemodel, as.character(i-1))$iter
		retval[i] = gtkTreeModelGetValue(treemodel,iter,as.integer(n))$value
	}
	
	return(retval)
}


fileChoose <- function(action="cat", text = "Select a file...", type="open", ...) 
{
	gfile(text=text, type=type, ..., action = action, handler =
		function(h,...) { do.call(h$action, list(h$file)) }, toolkit=guiToolkit("RGtk2")
		)
			
}


clearComboModel <- function(combo)
{
	gtkComboBoxSetActive(combo,-1)
	model = gtkComboBoxGetModel(StateEnv$itersCombo)
	Nelements = gtkTreeModelIterNChildren(model)
	for(i in 1:Nelements)
	{
		gtkComboBoxRemoveText(combo, 0)
	}
}


# Next four functions taken from rattle
Rtxt <- function(...)
{
  # Currently, on Windows we are waiting for 2.12.17 of  RGtk2 with
  # rgtk2_bindtextdomain().

#  if (.Platform$OS.type == "windows")
#    paste(...)
#  else
    gettext(paste(...), domain="R-WMCapacity")
}

# This is used to avoid the string being identified as a translation, as in
# RtxtNT(paste(vals ...))

RtxtNT <- Rtxt


packageIsAvailable <- function(pkg, msg=NULL)
{
  if (pkg %notin% rownames(installed.packages()))
  {
    if (not.null(msg))
      if (questionDialog(sprintf(Rtxt("The package '%s' is required to %s.",
                                      "It does not appear to be installed.",
                                      "A package can be installed",
                                      "with the following R command:",
                                      "\n\ninstall.packages('%s')",
                                      "\n\nThis will allow access to use",
                                      "the full functionality of %s.",
                                      "\n\nWould you like the package to be installed now?"),
                                 pkg, msg, pkg, crv$appname)))
      {
        install.packages(pkg)
        return(TRUE)
      }
    return(FALSE)
  }
  else
    return(TRUE)
}


"%notin%" <- function(x,y) ! x %in% y


.openCSVFile <- function(filename)
{
	dataset = try(read.csv(filename), silent=T)
	if (inherits(dataset, "try-error")){
		gmessage(paste("The file '",filename,"' could not be loaded. The probable cause this of error is that the file is not in CSV format.",
						sep=""), title="File error",icon="error",toolkit=guiToolkit("RGtk2"))
		return(0)
	}
	
	
	#.setDataForColumnSelection(dataset,filename=filename)
	
}

.loadProjectFile <- function(projectFile)
{
	myNewEnv = new.env(parent=globalenv())
	loaded = try(load(projectFile,env=myNewEnv), silent=T)
	if (inherits(loaded, "try-error")){
		gmessage(paste("The file '",projectFile,"' could not be loaded. The probable cause of this error is that the file is not an R save file.",
						sep=""), title="File error",icon="error",toolkit=guiToolkit("RGtk2"))
		return(0)
	}
	
	#.loadProject(myNewEnv$wommbatAnalysis,projectFile)	
}


.setInitialSensitivity<-function()
{
	.activeColumnSelection(FALSE)
	#.womActiveModelTab(FALSE)
	#.womActiveAnalysisTab(FALSE)
	#.womActiveDiagnosticsTab(FALSE)
	#.womActiveResultsTab(FALSE)
	#.womActiveSaveTab(FALSE)

}

