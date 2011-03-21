.columnsTreeCols <- function(name)
{
COLUMNS <- c(columnname=0,visible=1)
as.integer(COLUMNS[name])
}

.columnsOfInterestTreeCols <- function(name)
{
COLUMNS <- c(columnname=0,type=1,visible=2)
as.integer(COLUMNS[name])
}

.clicked_data_fixed_iv <- function(button)
{
	.addColumnOfInterest("fixed")
}

.clicked_data_random_iv <- function(button)
{
	.addColumnOfInterest("random")
}

.clicked_data_dv <- function(button)
{
	treeview = theWidget("treeviewDataColumns")
	model <- gtkTreeViewGetModel(treeview)
	selection = .getColumnNameSelection()
	
	# If nothing is selected, or a nonvisible row is selected, do nothing
	if(is.null(selection$iter) | selection$vis==FALSE){
		return()
	}
	
	# Get the name of the column
	selCol <- model$get(selection$iter, .columnsTreeCols("columnname"))[[1]]

	
	# Check column to see if it is numeric
	colnum = match(selCol,colnames(bfEnv$data))
		if(!is.na(colnum)){
			if(!is.numeric(bfEnv$data[,colnum]))
			{
				# Inappropriate data in column
				gmessage(paste("Invalid dependent variable data (must be numeric):",selCol), title="Data error",
					icon="error",toolkit=guiToolkit("RGtk2"))
				return(0)
			}
		}else{
			gmessage("Could not find selected column in data set! This should not have happened.",
				title="Interface error",icon="error",toolkit=guiToolkit("RGtk2"))
			return(0)
		}

	
	# If we had something there before, make it visible again
	oldText =  theWidget("entryDataY")$getText()
	if(oldText!=""){
			model$set(StateEnv$ColNamesIters[[oldText]], .columnsTreeCols("visible"), TRUE)
	}

	# Set the column name in the entry box
	theWidget("entryDataY")$setText(selCol)
		
	# Set the selected iter to be invisible
	model$set(selection$iter, .columnsTreeCols("visible"), FALSE)
	StateEnv$win$present()
	
}

.clicked_data_continuous_iv <- function(button)
{
	.addColumnOfInterest("continuous")
}

.clicked_data_remove_iv <- function(button)
{
	treeview1 = theWidget("treeviewDataColumns")
	model1 <- gtkTreeViewGetModel(treeview1)
	
	treeview2 = theWidget("treeviewDataIVs")
	model2 <- gtkTreeViewGetModel(treeview2)
	
	selection = .getColumnOfInterestSelection()
	
	# If nothing is selected, do nothing
	if(is.null(selection$iter)){
		return()
	}
	
	# Get the name of the column, make visible again
	selCol <- model2$get(selection$iter, .columnsOfInterestTreeCols("columnname"))[[1]]
	iter1 <- StateEnv$ColNamesIters[[selCol]]
	model1$set(iter1, .columnsTreeCols("visible"), TRUE)
	
	# Remove from treeview2
	model2$remove(selection$iter)
	
	# Set the selected iter to be invisible
	StateEnv$win$present()
}

.clicked_data_load_csv <- function(button)
{
	#freezeGUI()
	#on.exit(thawGUI())
	fileChoose(action=".OpenCSVFile", 
		type="open", 
		text="Select a CSV file...", 
		filter = list("CSV files" = list(patterns = c("*.csv","*.CSV")),
					  "text files" = list(mime.types = c("text/plain")),
					  "All files" = list(patterns = c("*")) 
					 )
		)
	
	
	#filename <- try(file.choose(), silent=T)
	#StateEnv$win$present()
	#if (inherits(filename, "try-error")) return()
	#.womOpenCSVFile(filename)
}


.OpenCSVFile <- function(filename)
{
	dataset = try(read.csv(filename), silent=T)
	if (inherits(dataset, "try-error")){
		gmessage(paste("The file '",filename,"' could not be loaded. The probable cause this of error is that the file is not in CSV format.",
						sep=""), title="File error",icon="error",toolkit=guiToolkit("RGtk2"))
		return(0)
	}
	
	.setDataForColumnSelection(dataset,filename=filename)
	
}


.clicked_data_load_dataframe <- function(button)
{
	freezeGUI()
	on.exit(thawGUI())
	.selectDataFrame(envir = .GlobalEnv)

}

.clicked_data_load_project <- function(button)
{
	#stuff here
}

.clicked_data_unlock <- function(button)
{
	#stuff here
}

.clicked_data_do_analysis <- function(button)
{
	treeview2 = theWidget("treeviewDataIVs")
	model2 <- gtkTreeViewGetModel(treeview2)

	# Check to make sure we defined at least one column of interest
	nRows <- gtkTreeModelIterNChildren(model2, NULL)	
	if(nRows==0)
	{
		gmessage(paste("You must define at least one column of interest."), title="Data error",
					icon="error",toolkit=guiToolkit("RGtk2"))
		return(0)
	}
	
	# Check to make sure we've defined necessary columns
	yText <- bfEnv$responseCol <- theWidget("entryDataY")$getText()
	
	if(yText=="")
	{
		gmessage(paste("You must define the dependent variable."), title="Data error",
					icon="error",toolkit=guiToolkit("RGtk2"))
		return(0)
	}

	# Disable the column selection
	.activeColumnSelection(FALSE)
	
	# Do the analysis here
	.buildAnalysisEnvironment()
	.doBayesFactorAnalysis()
	
	# Turn the page
	gtkNotebookSetCurrentPage(theWidget("notebook1"), .notebookPages("analysis"))

	
	StateEnv$win$present()
}

.buildAnalysisEnvironment <- function()
{
	analysisTypeCombo = StateEnv$typeCombo
	
	bfEnv$settings = list()
	bfEnv$settings$iterations = as.numeric(theWidget("entrySettingsIterations")$getText())
	bfEnv$settings$analysisType = analysisTypeCombo$getActiveText()
	
	model <- gtkTreeViewGetModel(theWidget("treeviewDataIVs"))
	dataset <- bfEnv$data
	
	# get the selected columns and types
	bfEnv$selectedCols <- .treeModelGetNthCol(model,.columnsOfInterestTreeCols("columnname"))
	bfEnv$selectedColsType <- .treeModelGetNthCol(model,.columnsOfInterestTreeCols("type"))	

	bfEnv$selectedColsType[bfEnv$selectedColsType=="continuous"] <- "C"
	bfEnv$selectedColsType[bfEnv$selectedColsType=="fixed"] <- "F"
	bfEnv$selectedColsType[bfEnv$selectedColsType=="random"] <- "R"
	
	.shapeDataForAnalysis()
	
	bfEnv$nFixedFac = nFac = sum(bfEnv$selectedColsType=="F")
	bfEnv$designMatrices = list()
	bfEnv$designMatrices[[2^nFac]] = matrix(nrow=0,ncol=0)	
	bfEnv$totalN = length(bfEnv$y)

}

.shapeDataForAnalysis <- function()
{
	data = bfEnv$data
	cols = colnames(data)
	selCols = bfEnv$selectedCols
	selType = bfEnv$selectedColsType
	
	bfEnv$y = data[,cols==bfEnv$responseCol]
	
	dataRandom = data.frame(row.names=1:length(bfEnv$y))
	dataFixed = data.frame(row.names=1:length(bfEnv$y))
	dataContinuous = data.frame(row.names=1:length(bfEnv$y))

	
	for(i in 1:length(selCols))
	{
		name = selCols[i]
		col = data[,cols==name]
		
		if(selType[i]=="F")
		{
			dataFixed = data.frame(dataFixed,as.factor(col)[drop=TRUE])
		}
		if(selType[i]=="R")
		{
			dataRandom = data.frame(dataRandom,as.factor(col)[drop=TRUE])
		}
		
		if(selType[i]=="C")
		{
			dataContinuous = data.frame(dataContinuous,as.factor(col)[drop=TRUE])
		}		
	}
	
	if(any(selType=="R")){
		colnames(dataRandom) = selCols[selType=="R"]
		bfEnv$dataRandom = dataRandom
	}else{
		bfEnv$dataRandom = NULL
	}
	if(any(selType=="F")){
		colnames(dataFixed) = selCols[selType=="F"]
		bfEnv$dataFixed = dataFixed
	}else{
		bfEnv$dataFixed = NULL
	}
	if(any(selType=="C")){
		colnames(dataContinuous) = selCols[selType=="C"]
		bfEnv$dataContinuous = dataContinuous
	}else{
		bfEnv$dataContinuous = NULL
	}
}

.doBayesFactorAnalysis <- function()
{
	analysisType = bfEnv$settings$analysisType
	
	if(analysisType=="All"){
		bfEnv$bayesFactors = 
				c(null=0,
				all.Nways(bfEnv,iterations=bfEnv$settings$iterations)
			   )*log10(exp(1))
	}
	
	.setAnalysisResults()
}


.changed_data_analysis_type <- function(button)
{
	#stuff here
}

.addColumnOfInterest <- function(type)
{
	treeview1 = theWidget("treeviewDataColumns")
	model1 <- gtkTreeViewGetModel(treeview1)
	
	treeview2 = theWidget("treeviewDataIVs")
	model2 <- gtkTreeViewGetModel(treeview2)
	
	selection = .getColumnNameSelection()
	
	# If nothing is selected, or a nonvisible row is selected, do nothing
	if(is.null(selection$iter) | selection$vis==FALSE){
		return()
	}
	
	# Get the name of the column
	selCol <- model1$get(selection$iter, .columnsTreeCols("columnname"))[[1]]
	
	# Add the column name in treeview2
	iter <- model2$append(NULL)$iter
	model2$set(iter, 
			  .columnsOfInterestTreeCols("columnname"), selCol,
			  .columnsOfInterestTreeCols("type"), type,
			  .columnsOfInterestTreeCols("visible"), TRUE)
		
	
		
	# Set the selected iter to be invisible
	model1$set(selection$iter, .columnsTreeCols("visible"), FALSE)
	StateEnv$win$present()
}

.getColumnNameSelection <- function()
{
	treeview = theWidget("treeviewDataColumns")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	

	vis <- model$get(iter, .columnsTreeCols("visible"))[[1]]
	return(list(iter=iter,visible=vis))
}

.setDataForColumnSelection <- function(dataset,removeModels=TRUE,filename="")
{
	#cat("Setting columns\n")
	
	# Confirm that it is ok to reset the interface.
	if(length(bfEnv$designMatrices)>0)
	{
		confirm = gconfirm("This will remove all existing models and analyses.", title="Reset interface?",toolkit=guiToolkit("RGtk2"))
		if(!confirm) return(0)
	}
	
	# check to make sure there are column names
	cols <- colnames(dataset)
	if(any(is.null(cols)))
	{
		gmessage("The data set could not be loaded, because there were no column names.", title="Data error",icon="error",toolkit=guiToolkit("RGtk2"))
		return(0)
	}
	
	# check to make sure that there are enough columns and rows.
	myDims = dim(dataset)
	if(myDims[1]<2 | myDims[2]<2)
	{
		gmessage("The data set could not be loaded, because too few rows (<2) or columns (<2).", title="Data error",icon="error",toolkit=guiToolkit("RGtk2"))
		return(0)

	}
	
	# check to make sure all column names are unique. If not, error.
	if(any(table(cols))>1){
	gmessage("The data set could not be loaded, because some column names were duplicates.", title="Data error",icon="error",toolkit=guiToolkit("RGtk2"))
		return(0)
	}
	bfEnv$data <- dataset
	theWidget("entryDataFilename")$setText(filename)
	
	#cat("creating treeviews\n")
	
	.createNewColumnsOfInterestTreeModel()
	.createNewColumnSelectionTreeModel()
	.createColumnSelectionColumns()
	.createColumnsOfInterestColumns()
	
	treeview1 <- theWidget("treeviewDataColumns")
	model1 <- gtkTreeViewGetModel(treeview1)
	
	treeview2 <- theWidget("treeviewDataIVs")
	model2 <- gtkTreeViewGetModel(treeview2)
	

	StateEnv$ColNamesIters=list()
	
	# Add dataset column names to treeview model
	for(col in cols){
		iter <- model1$append(NULL)$iter
		model1$set(iter,
	  		  .columnsTreeCols("columnname"), col,
			  .columnsTreeCols("visible"), TRUE)
		StateEnv$ColNamesIters[[col]] = iter
	}
	
	theWidget("entryDataY")$setText("")
	
	.activeColumnSelection(TRUE)
	.activeSaveBox(TRUE)
	.clearInterfaceForDataLoad(removeModels)
	#StateEnv$win$present()	
}

.clearInterfaceForDataLoad <- function(removeModels=TRUE)
{
	# Models
	if(removeModels){
		bfEnv$designMatrices=NULL
		#bfEnv$data = NULL
		bfEnv$dataRandom = NULL
		bfEnv$dataFixed = NULL
		bfEnv$selectedCols = NULL
		bfEnv$selectedColsType = NULL
		bfEnv$bayesFactors = NULL
		bfEnv$y = NULL
		bfEnv$totalN = NULL
		bfEnv$nFixedFac = NULL
		bfEnv$responseCol = NULL
	}
	
	# Analysis
	if(!is.null(theWidget("treeviewAnalysisBFs")$getModel()))
		{
			theWidget("treeviewAnalysisBFs")$getModel()$clear()
		}		
	.activeAnalysisTab(FALSE)	
	
}


.activeColumnSelection<-function(status=TRUE)
{

	theWidget("boxDataColumns")$setSensitive(status)
	theWidget("boxDataDoAnalysis")$setSensitive(status)

}

.createNewColumnsOfInterestTreeModel <- function()
{
	model2 <- gtkTreeStoreNew("gchararray","gchararray","gboolean")
	treeview2 <- theWidget("treeviewDataIVs")
	
	gtkTreeViewSetModel(treeview2, model2)
}


.createNewColumnSelectionTreeModel <- function()
{
	model1 <- gtkTreeStoreNew("gchararray","gboolean")
	treeview1 <- theWidget("treeviewDataColumns")
	
	gtkTreeViewSetModel(treeview1, model1)
}

.createColumnsOfInterestColumns <- function()
{
	treeview2 <- theWidget("treeviewDataIVs")
	model2 <- gtkTreeViewGetModel(treeview2)
	


	#Create right columns tree view
	if(is.null(gtkTreeViewGetColumn(treeview2,0))){
		
		# column name
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.0)

		col.offset <- treeview2$insertColumnWithAttributes(-1, "Column", renderer, 
  								text = .columnsOfInterestTreeCols("columnname"),
								visible = .columnsOfInterestTreeCols("visible"))
								
		column <- treeview2$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# column type
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.0)

		col.offset <- treeview2$insertColumnWithAttributes(-1, "Type", renderer, 
  								text = .columnsOfInterestTreeCols("type"),
								visible = .columnsOfInterestTreeCols("visible"))
								
		column <- treeview2$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
	}

}

.createColumnSelectionColumns <- function()
{
	treeview1 <- theWidget("treeviewDataColumns")
	model1 <- gtkTreeViewGetModel(treeview1)
	
	#Create left columns tree view
	if(is.null(gtkTreeViewGetColumn(treeview1,0))){
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.0)

		col.offset <- treeview1$insertColumnWithAttributes(-1, "Column name", renderer, 
  								text = .columnsTreeCols("columnname") ,
								visible = .columnsTreeCols("visible"))
								
		column <- treeview1$getColumn(col.offset - 1)
		column$setClickable(TRUE)
	}
}


.getColumnOfInterestSelection<-function()
{
	treeview = theWidget("treeviewDataIVs")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	

	return(list(iter=iter))
}


.setupDataAnalysisTypeComboBox<-function()
{
	typeSpace = theWidget('boxDataAnalysisType')
	StateEnv$typeCombo = gtkComboBoxNewText()
	gtkComboBoxAppendText(StateEnv$typeCombo, "All")
	typeSpace$packStart(StateEnv$typeCombo,FALSE,FALSE,0)
	gtkComboBoxSetActive(StateEnv$typeCombo,0)
	# Connect parameter type signal
	StateEnv$handlers$dataAnalysisTypeCombo <- gSignalConnect(StateEnv$typeCombo, "changed", .changed_data_analysis_type)
}
