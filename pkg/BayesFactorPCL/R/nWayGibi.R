

AOVnWayGibi <- function() {
  list(


    get.variables = function() {
      c("Dependent variable","Participants")
    },


    get.columns = function() {
      c("Independent variable")
    },


    assign.variables = function(dataframe, vars, cols, result) {

      #result$ready <- ( ! ("" %in% unlist(vars))) && (length(unlist(cols)) > 0)
      result$ready <- ( ! ("" %in% vars[["Dependent variable"]])) && (length(unlist(cols)) > 0)

      result
    },


    create.effects = function(vars, cols) {

      factors <- unname(unlist(cols))
      effects <- list("Available"=as.list(factors))
	  
	  
	  
#      if (length(factors) > 1) {
#        for (i in 2:length(factors)) {
#          nway <- c()
#          combinations <- combn(factors, i)
#          for (j in 1:dim(combinations)[2]) {
#            desc <- paste(combinations[,j], collapse=" * ")
#            nway <- append(nway, desc)
#          }
#          effects[[paste(i, "way")]] <- as.list(nway)
#        }
#      }



      effects
    },


    create.model = function() {
      GIBI::Model$new(

        level.one = c("Factors"),
		
		# don't need level two
        level.two = NULL,#c("DON'T NEED THIS"),

        priors    = create.priors(),
        analysis.params = create.analysis.params()
      )
    },

    create.models.table.headings = function() {

      c("#factors", "#models", "delta scale")
      
    },

    create.models.table.row = function(model) {


      level.one <- model$get.level.one()
      priors    <- model$get.priors()

	  nfactors = length(unlist(level.one[["Factors"]]))	
      nmodels = 2^(2^nfactors-1)
	  
      delta.scale <- priors$get("delta (effect size)")$get("Cauchy scale")$get.value()
 

      c(nfactors, nmodels, delta.scale)

    },

    create.analysis.table.headings = function() {
      c("Iterations","Subsets","#analyses")
    },

    create.analysis.table.row = function(model) {

      params <- model$get.analysis.params()
	  
	  iterations = params$get("(MC)MC")$get("Iterations")$get.value()
      subsets <- params$get("General ANOVA")$get("Model subsets")$get.value()
  	  
  	  nfactors <- length(unlist(model$get.level.one()[["Factors"]]))
	  if(subsets=="All"){
      	nmodels <- 2^(2^nfactors-1)
      }else{
      	nmodels <- 2^nfactors
      }

      c(iterations,subsets,nmodels)

      },

    perform.analysis = function(dataframe, variables, model, progress.callback) {
	  
	  priors    <- model$get.priors()
      delta.scale <- as.numeric(priors$get("delta (effect size)")$get("Cauchy scale")$get.value())

	  params <- model$get.analysis.params()
      iterations = params$get("(MC)MC")$get("Iterations")$get.value()
      burnin =  params$get("(MC)MC")$get("Burnin Iterations")$get.value()
      subsets <- params$get("General ANOVA")$get("Model subsets")$get.value()
      
      if(subsets=="All"){
      	only.top=FALSE
      }else{
      	only.top=TRUE
      }
      ivs <- unlist(model$get.level.one()[["Factors"]])
      
      dv <- dataframe[,variables[["Dependent variable"]]]
      fixed <- dataframe[,ivs]
	  if(variables[["Participants"]]==""){
	  	subs <-NULL
	  }else{
	    subs <- dataframe[,variables[["Participants"]]]
      }
      
      
 	  results<-all.Nways(dv,dataFixed=fixed,dataRandom=subs,iterations = iterations,samples=TRUE,only.top=only.top,rscale=delta.scale)
	  progress.callback(model, 100)

      return(results)
    },

    get.diagnostics.slice.count = function(model) {
      length(model$get.results()[[1]])
    },

    get.diagnostics.graph.count = function() {
      1
    },

    create.diagnostics.graph = function(model, n, slice) {
      
      chains <- model$get.results()[[2]][[slice]]
      modname <- names(model$get.results()[[1]])[slice]
      
	  if(length(chains)==1){
	  	plot(1,1,typ='n', main="null")
	  	text(1,1,"NA")
	  }else{
	    cummeans <- logCumMeanExpLogs(chains)
	  	
	  	plot(cummeans,ty='l',main=modname, ylab="Cum. Mean",xlab="Iteration")
	  	
	  	#plot(log10(1:length(cummeans)), logCumMeanExpLogs(chains),ty='l',main=modname,ylab="Cum. Mean",xlab="Iteration",axes=FALSE)
	  	#axis(1,at=0:floor(log10(length(cummeans))),10^(0:floor(log10(length(cummeans)))))
	  	#box()
	  }

    },

    create.results.table.headings = function() {
      c("#models")
    },

    create.results.table.row = function(model) {

      results <- model$get.results()

      c(length(results[[1]]))

    },

	create.parameters.table.headings = function() {

      c("Model", "Bayes factor", "log10 BF")

	},

	create.parameters.table.rows = function(model) {
	
	  res <- model$get.results()[[1]]
      ests <- cbind(names(res), 10^res,res)
      dimnames(ests) <- NULL
      listEsts <- apply(ests,1,as.list)
	  
	  return(listEsts)
	}

  )
}
