BFInfo <- function()
{
	print("Package BayesFactor")
	print(packageDescription("BayesFactor")$Version)
  myRev <- '$Rev$'
  myDate <- '$Date$'
  nRev <- nchar(myRev)
  nDate <- nchar(myDate)
  myRev <- substr(myRev,2,nRev-1)
  myDate <- substr(myDate,2,nDate-1)
  print(myRev)
	print(myDate)
  invisible()
}