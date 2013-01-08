setGeneric("%same%", function(x, y) standardGeneric("%same%"))

setGeneric("compare", function(numerator, denominator, data, ...) standardGeneric("compare"))

setGeneric("recompute", function(x, progress=FALSE, ...) standardGeneric("recompute"))

setGeneric("posterior", function(model, index, data, iterations, ...) standardGeneric("posterior"))

setGeneric("extractBF", function(x, ...) standardGeneric("extractBF"))
