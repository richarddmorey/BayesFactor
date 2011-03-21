
.activeAnalysisTab<-function(status=TRUE)
{

	theWidget("boxAnalysisPage")$setSensitive(status)

}

.analysisTreeCols <- function(name)
{
COLUMNS <- c(number=0,name=1,df=2,terms=3,numlogbf=4,logbf=5,bf=6,iterations=7,mcmcse=8,analyze=9,plot=10)
as.integer(COLUMNS[name])
}


.setAnalysisResults <- function()
{
	bfs = bfEnv$bayesFactors
	
	bf.df = data.frame(modelNum=1:length(bfs)-1,names=names(bfs),log10bf=bfs,iters=bfEnv$settings$iterations)
	
	.createAnalysisTreeview()
	#treeview <- theWidget("treeviewAnalysisBFs")
	#treemodel <- gtkTreeViewGetModel(treeview)
	treemodel <- StateEnv$analysisModel
	
	for(i in 1:length(bfs))
	{
		iter <- treemodel$append(NULL)$iter
		treemodel$set(iter, 
			  .analysisTreeCols("number"), bf.df[i,1],
			  .analysisTreeCols("name"), bf.df[i,2],
			  .analysisTreeCols("df"), 0,
			  .analysisTreeCols("terms"), 0,
			  .analysisTreeCols("numlogbf"), bf.df[i,3],
			  .analysisTreeCols("logbf"), .saneNum(bf.df[i,3],2),
			  .analysisTreeCols("bf"), .logbfTobfText(bf.df[i,3]),
			  .analysisTreeCols("iterations"), bf.df[i,4],
			  .analysisTreeCols("mcmcse"), 0,
			  .analysisTreeCols("analyze"), FALSE,
			  .analysisTreeCols("plot"), FALSE
			  )

	}
	
	bfEnv$bfDataFrame = bf.df
	
	StateEnv$win$present()

}

.createAnalysisTreeStore <- function()
{
	model <- gtkTreeStoreNew("gint","gchararray","gint","gint","gdouble","gchararray","gchararray","gint","glong","gboolean","gboolean")
	treeview <- theWidget("treeviewAnalysisBFs")
	
	StateEnv$analysisModel = model
	
	analysisSortModel <- gtkTreeModelSortNewWithModel(child.model = model)
	gtkTreeViewSetModel(treeview, analysisSortModel)
	
}


.createAnalysisTreeview <- function()
{
	treeview <- theWidget("treeviewAnalysisBFs")
	
	#Create tree view
	if(is.null(gtkTreeViewGetColumn(treeview,0))){
			
		# number
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Number", renderer, 
  								text = .analysisTreeCols("number"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.analysisTreeCols("number"))
		
		# name
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Name", renderer, 
  								text = .analysisTreeCols("name"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# log Bayes factor
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "LogBF", renderer, 
  								text = .analysisTreeCols("logbf"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.analysisTreeCols("numlogbf"))
		
		# Bayes factor
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "BF", renderer, 
  								text = .analysisTreeCols("bf"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.analysisTreeCols("numlogbf"))
		
	}
}
