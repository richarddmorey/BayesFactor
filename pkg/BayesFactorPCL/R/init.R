
BFAnovaInit<-function()
{
  tclRequire("BWidget")
} 

BFAnovaFonts=function(){
list(
fontHeading = tkfont.create(family="times",size=16,weight="bold")
)
}

BFAnovaHelp=function(){
list(
ColSelHelp = "Empty",

EffSelHelp="Empty", 

PriorSetHelp="Empty",

MCMCHelp="Empty",


InstrHelp="Click help at any time for help regarding the options in a screen.",

OutputHelp="Empty",

InitInstructions="Instructions go here."
)			
}
