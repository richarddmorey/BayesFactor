
.activeAnalysisTab<-function(status=TRUE)
{

	theWidget("boxAnalysisPage")$setSensitive(status)

}

.setAnalysisResults <- function()
{
	bfs = bfEnv$bayesFactors
	
	bf.df = data.frame(modelNum=1:length(bfs)-1),names=names(bfs),log10bf=bfs,iters=bfEnv$settings$iterations)
	

}

.createAnalysisTreeStore <- function()
{
	model <- gtkTreeStoreNew("gchararray","gchararray","gchararray","gboolean")
	treeview <- theWidget("treeviewAnalysisBFs")
	
	analysisSortModel <- gtkTreeModelSortNewWithModel(child.model = model)
	gtkTreeViewSetModel(treeview, analysisSortModel)
	
	
}


.createAnalysisTreeview <- function()
{
	treeview <- theWidget("treeviewAnalysisBFs")
	
	#Create tree view
	if(is.null(gtkTreeViewGetColumn(treeview,0))){
		
		
		# name
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Name", renderer, 
  								text = .womDefinedModelsTreeCols("name"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# hasResults
		renderer <- gtkCellRendererToggleNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Results?", renderer, 
  								active = .womDefinedModelsTreeCols("hasResults")
								)
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(FALSE)
		
		
		# time
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Time", renderer, 
  								text = .womDefinedModelsTreeCols("timeAnalyzed"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)

				
		# acc rate
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Acc. Rate", renderer, 
  								text = .womDefinedModelsTreeCols("accRate"),
								foreground = .womDefinedModelsTreeCols("accColor"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.womDefinedModelsTreeCols("accRate"))
		
		
		# iterations
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Iterations", renderer, 
  								text = .womDefinedModelsTreeCols("iterations"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.womDefinedModelsTreeCols("iterations"))
		
		# effective iterations
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Effective Iterations", renderer, 
  								text = .womDefinedModelsTreeCols("effectiveIterations"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.womDefinedModelsTreeCols("effectiveIterations"))
		
		# burnin
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Burnin", renderer, 
  								text = .womDefinedModelsTreeCols("burnin"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.womDefinedModelsTreeCols("burnin"))
		
		# hybrid EPS 
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Lower eps.", renderer, 
  								text = .womDefinedModelsTreeCols("epsLow"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)

		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Upper eps.", renderer, 
  								text = .womDefinedModelsTreeCols("epsUpp"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# Leapfrog	
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Leapfrog", renderer, 
  								text = .womDefinedModelsTreeCols("leapfrog"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# use Metropolis
		renderer <- gtkCellRendererToggleNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "MH instead?", renderer, 
  								active = .womDefinedModelsTreeCols("useMet")
								)
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(FALSE)
	
		# metrop scale
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "MH scale", renderer, 
  								text = .womDefinedModelsTreeCols("metropScale"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# metrop scale
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "MH thin", renderer, 
  								text = .womDefinedModelsTreeCols("metropThin"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
	}
}
