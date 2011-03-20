# this stores application (non-project) state information
if (!exists("StateEnv", environment(), inherits=FALSE)) {
	StateEnv <- new.env()
}
# stores all project information
if (!exists("wommbatAnalysis", environment(), inherits=FALSE)) {
	bfEnv <- new.env()
}

BFgui <- function(project = NULL, projectFile= NULL,CSVfile = NULL, dataFrame = NULL) {

	options("guiToolkit"="RGtk2")

	StateEnv$update <- list()
	#bfEnv <- list()

	StateEnv$echo.to.console <- TRUE
	StateEnv$echo.to.log <- TRUE
	StateEnv$Graphics <- list()
	
	StateEnv$GUI <- gtkBuilder()
	filename <- getpackagefile("AOVinterface.ui")
	res <- StateEnv$GUI$addFromFile(filename)
	if (!is.null(res$error))
		stop("ERROR: ", res$error$message)
		
	StateEnv$win <- theWidget("window1")
	StateEnv$win$setTitle("Bayes Factors")

	# connect the callbacks (event handlers)
	#gladeXMLSignalAutoconnect(StateEnv$GUI)
	StateEnv$GUI$connectSignals(NULL)
	StateEnv$handlers = list()

	#.loadTooltips()
	#.populateTextviews()
	
	#.setInitialSensitivity()

	.setupDataAnalysisTypeComboBox()
	
		# load files/dataframes
	if(!is.null(dataFrame)){
		if(is.data.frame(dataFrame)){
			theWidget("entryDataFilename")$setText("<Loaded from dataframe>")
			.setDataForColumnSelection(dataFrame)
		}
	}else if(!is.null(CSVfile))
	{
		.OpenCSVFile(CSVfile)
	}else if(!is.null(project))
	{
		#.loadProject(project)
	}else if(!is.null(projectFile))
	{
		#.loadProjectFile(projectFile)
	}
	
	
	StateEnv$win$present()
	return(invisible(NULL))
}



