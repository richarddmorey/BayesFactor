
create.priors <- function() {
  Params$new(params=list(
    Params$new("delta (effect size)", list(
      Number$new("Cauchy scale", 1))
      )
    )
  )
}

create.analysis.params <- function() {
  Params$new(params=list(
  	Params$new("General ANOVA", list(
  	 List$new("Model subsets", c("Top", "All"), "Top"))
  	),
    Params$new("(MC)MC", list(
      Number$new("Iterations", 10000),
      Number$new("Burnin Iterations", 200))
	  )
	)
  )
}

AOV1Gibi <- function() {
  list(


    get.variables = function() {
      c("Response variable")
    },


    get.columns = function() {
      c("Group variable")
    },


    assign.variables = function(dataframe, vars, cols, result) {

      result$ready <- ( ! ("" %in% unlist(vars))) && (length(unlist(cols)) > 0)

      result
    },


    create.effects = function(vars, cols) {

      factors <- unname(unlist(cols))
      effects <- list("Groups"=as.list(factors))

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

        level.one = c("IVs"),
		
		# don't need level two
        level.two = NULL,#c("DON'T NEED THIS"),

        priors    = create.priors(),
        analysis.params = create.analysis.params()
      )
    },

    create.models.table.headings = function() {

      c("IV", "delta scale")
      
    },

    create.models.table.row = function(model) {


      level.one <- model$get.level.one()
      priors    <- model$get.priors()

      iv <- level.one[["IVs"]][1]

      delta.scale <- priors$get("delta (effect size)")$get("Cauchy scale")$get.value()
 

      c(iv, delta.scale)



    },

    create.analysis.table.headings = function() {
      c("Iterations", "Burnin")
    },

    create.analysis.table.row = function(model) {

      params <- model$get.analysis.params()

      c(params$get("(MC)MC")$get("Iterations")$get.value(),
        params$get("(MC)MC")$get("Burnin Iterations")$get.value()
        )

      },

    perform.analysis = function(dataframe, variables, model, progress.callback) {
	  
	  priors    <- model$get.priors()
      delta.scale <- as.numeric(priors$get("delta (effect size)")$get("Cauchy scale")$get.value())
 
	  
	  params <- model$get.analysis.params()
      
      iv <- model$get.level.one()[["IVs"]][[1]]
	  ivCol <- dataframe[,iv]
      
      dv <- dataframe[,variables[["Response variable"]]]
      
      # Convert into proper matrix form
      maxN = max(table(ivCol))
      unpadded <- tapply(dv,ivCol,c)
      padded <- lapply(unpadded,function(v,n) c(v,rep(NA,n-length(v))),n=maxN)
      datmat <- matrix(unlist(padded),nrow=maxN)
      myNames <- names(padded)
      
      iterations = params$get("(MC)MC")$get("Iterations")$get.value()
      burnin =  params$get("(MC)MC")$get("Burnin Iterations")$get.value()
      
      
      gibi.callback = function(percent){
      	progress.callback(model, percent)
      }
      results = oneWayAOV.Gibbs(y=datmat,
      							iterations=iterations,
      							rscale=delta.scale,
      							gibi=gibi.callback)
      

	  
	  # Make nice chain names
	  chainNames <- colnames(results[[1]])
	  chainNames[1+1:length(myNames)]<-myNames
	  colnames(results[[1]])<-chainNames
  
      
      return(list(results))

    },

    get.diagnostics.slice.count = function(model) {
      #ncol(model$get.results()[[1]][[1]])
      ncol(globalenv()$resultsList[[1]][[1]]) 
    },

    get.diagnostics.graph.count = function() {
      3
    },

    create.diagnostics.graph = function(model, n, slice) {
      

      
      #chains <- model$get.results()[[1]][[1]]
       chains <- globalenv()$resultsList[[1]][[1]]
  
  	   parname <- colnames(chains)[slice]
      
      if(n==1){
      	plot(as.vector(chains[,slice]), type="l", main = paste(parname,"chain"), ylab="parname",xlab="Iteration")
      }else if(n==2){
      	plot(density(as.vector(chains[,slice])), main = paste(parname,"kernel density"),xlab="parname",ylab="")
      }else if(n==3){
      	plot(1,1)
      	plot(acf(as.vector(chains[,slice])), main = paste(parname," ACF"))
      }
    },

    create.results.table.headings = function() {
      c("Bayes factor","log Bayes factor")
    },

    create.results.table.row = function(model) {

	  #results <- model$get.results()	
      results <- globalenv()$resultsList
      results[[1]][[3]]<-log(results[[1]][[2]])
      c(results[[1]][[2]],results[[1]][[3]])

    },

	create.parameters.table.headings = function() {

      c("Parameter", "Mean")

	},

	create.parameters.table.rows = function(model) {
	  
	  #chains <- model$get.results()[[1]][[1]]
	  chains <- globalenv()$resultsList[[1]][[1]]
      ests <- cbind(colnames(chains), colMeans(chains))
      rownames(ests) <- NULL
      listEsts <- apply(ests,1,as.list)
	  return(listEsts)
	}

  )
}
