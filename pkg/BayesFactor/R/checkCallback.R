

checkCallback = function(callback, ... ){
  ret = callback(...)
  if(ret) stop("Operation cancelled by callback function. ", ret)
  return(ret)
}