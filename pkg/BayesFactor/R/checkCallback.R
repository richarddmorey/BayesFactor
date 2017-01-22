

checkCallback = function(callback, ... ){
  ret = as.integer(callback(...))
  if(ret) stop("Operation cancelled by callback function. ", ret)
  return(ret)
}
