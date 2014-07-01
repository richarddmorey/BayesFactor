#' Plot a sequence of Bayes factors
#' 
#' This function creates a lineplot of a sequence of log Bayes factors. As a default, the function expects raw Bayes factors (i.e., non-logged Bayes factors). #' If you provide Bayes factors that are already logged (e.g., from the output of the ttest.tstat function), set \code{log.it} to \code{FALSE}.
#' This function is in particular useful for plotting the trajectory of a sequential Bayes factor test
#' @title Plot a Bayes factor object
#' @param x A vector o f numbers for the x axis
#' @param BF A vector of Bayes factors (same length as x)
#' @param xlab Label for x axis
#' @param ylab Label for y axis
#' @param log.it Should the log of the Bayes factor be used for plotting?
#' @param forH1 If \code{TRUE}, positive BFs mean evidence in favor of H1. This is the default in the Bayes factor package.
#' @author Felix Schoenbrodt (\email{felix@nicebread.de})
#' @examples
#' \dontrun{
#' ## Sleep data from t test example
#' data(sleep)
#' 
#' # Compute accumulating evidence from n1=5 participants to n2=10 participants
#' BF <- c()
#' for (i in 5:10) {
#' 	BF0 <- ttestBF(x = sleep$extra[sleep$group==1][1:i], y = sleep$extra[sleep$group==2][1:i], paired=TRUE)
#' 	BF <- c(BF, as.vector(BF0))
#' }
#' 
#' BFplot(5:10, BF)
#' }



seqBFplot <- function(x, BF, xlab="n", ylab="log(BF)", log.it=TRUE, forH1=TRUE) {
	if (length(x) != length(BF)) stop("`x` and `BF` should have the same length")
	if (length(x) == 1) stop("`x`and `BF` must habe length > 1")
		
	if (log.it==TRUE) BF <- log(BF)
		
	df <- data.frame(x, BF)
	p1 <- ggplot(df, aes(x=x, y=BF)) + geom_line()+ theme_bw() + ylab(ylab) + xlab(xlab)
	p1 <- p1 + geom_hline(yintercept=c(c(-log(c(100, 30, 10, 3)), log(c(3, 10, 30, 100)))), linetype="dotted", color="darkgrey")
	p1 <- p1 + geom_hline(yintercept=log(1), linetype="dashed", color="darkgreen")

	p1 <- p1 + annotate("text", x=max(x)*1.5, y=-5.15, label=paste0("Extreme~H[", ifelse(forH1==TRUE,0,1), "]"), 
		hjust=1, vjust=.5, size=3, color="black", parse=TRUE)
	p1 <- p1 + annotate("text", x=max(x)*1.5, y=-4.00, label=paste0("Very~strong~H[", ifelse(forH1==TRUE,0,1), "]"), 
		hjust=1, vjust=.5, size=3, color="black", parse=TRUE)
	p1 <- p1 + annotate("text", x=max(x)*1.5, y=-2.85, label=paste0("Strong~H[", ifelse(forH1==TRUE,0,1), "]"), 
		hjust=1, vjust=.5, size=3, color="black", parse=TRUE)
	p1 <- p1 + annotate("text", x=max(x)*1.5, y=-1.7 , label=paste0("Moderate~H[", ifelse(forH1==TRUE,0,1), "]"), 
		hjust=1, vjust=.5, size=3, color="black", parse=TRUE)
	p1 <- p1 + annotate("text", x=max(x)*1.5, y=-.55 , label=paste0("Anectodal~H[", ifelse(forH1==TRUE,0,1), "]"), 
		hjust=1, vjust=.5, size=3, color="black", parse=TRUE)

	p1 <- p1 + annotate("text", x=max(x)*1.5, y=5.15, label=paste0("Extreme~H[", ifelse(forH1==TRUE,1,0), "]"), 
		hjust=1, vjust=.5, size=3, color="black", parse=TRUE)
	p1 <- p1 + annotate("text", x=max(x)*1.5, y=4.00, label=paste0("Very~strong~H[", ifelse(forH1==TRUE,1,0), "]"), 
		hjust=1, vjust=.5, size=3, color="black", parse=TRUE)
	p1 <- p1 + annotate("text", x=max(x)*1.5, y=2.86 , label=paste0("Strong~H[", ifelse(forH1==TRUE,1,0), "]"), 
		hjust=1, vjust=.5, size=3, color="black", parse=TRUE)
	p1 <- p1 + annotate("text", x=max(x)*1.5, y=1.7  , label=paste0("Moderate~H[", ifelse(forH1==TRUE,1,0), "]"), 
		hjust=1, vjust=.5, size=3, color="black", parse=TRUE)
	p1 <- p1 + annotate("text", x=max(x)*1.5, y=.55  , label=paste0("Anectodal~H[", ifelse(forH1==TRUE,1,0), "]"), 
		hjust=1, vjust=.5, vjust=.5, size=3, color="black", parse=TRUE)

	# set scale ticks
	p1 <- p1 + scale_y_continuous(breaks=c(c(-log(c(100, 30, 10, 3)), 0, log(c(3, 10, 30, 100)))), labels=c("-log(100)", "-log(30)", "-log(10)", "-log(3)", "log(1)", "log(3)", "log(10)", "log(30)", "log(100)"))
	
	# custom labeler: find breaks with pretty nbumers, and not more than 5
	myBreaks <- function(y){
		# access variable from parent environment
		MIN <- min(get("x", parent.frame()))
		MAX <- max(get("x", parent.frame()))

		# find good divisor
		steps <- c(2, 4, 5, 10, 15, 20, seq(30, 10000, by=10))
		i <- 1
		repeat {
			mod <- (MAX-MIN+1) %/% steps[i]
			if (mod <= 5) {break} else {i <- i+1}
		}
		# x is taken from the parent.env
	    breaks <- seq(steps[i], MAX, by=steps[i])
	    names(breaks) <- breaks
	    breaks
	}
	p1 <- p1 + scale_x_continuous(breaks=myBreaks)
	return(p1)
}