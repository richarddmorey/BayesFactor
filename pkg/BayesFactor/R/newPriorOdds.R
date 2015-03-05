#'Create prior odds from a Bayes factor object
#'
#'Create a prior odds object from a Bayes factor object
#'
#'This function takes a Bayes factor object and, using its structure and
#'specified type of prior odds, will create a prior odds object.
#'
#'For now, the only type is "equal", which assigns equal prior odds to all 
#'models.
#'
#'@param bf A BFBayesFactor object, eg, from an analysis
#'@param type The type of prior odds to create (by default "equal"; see details)
#'@return A (prior) BFodds object, which can then be multiplied by the 
#'BFBayesFactor object to obtain posterior odds.
#'@author Richard D. Morey (\email{richarddmorey@@gmail.com})
#'@keywords misc
#'@export
newPriorOdds = function(bf, type = "equal"){
  BFodds(bf)
}