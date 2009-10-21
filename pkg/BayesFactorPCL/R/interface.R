
BFAnovaInitInstructions <- function(){	
	thisEnv=environment()
	
	tt<-tktoplevel()
	done <- tclVar(0) # Are we done yet?
	tkgrab.set(tt)
	tktitle(tt) = "Instructions"
	frameOuter <- tkframe(tt)
	tkgrid(tklabel(frameOuter,text="BF-ANOVA Instructions",font=BFAnovaFonts()$fontHeading),row=0,columnspan=3)
	
	# Instructions frame
	frameInstructions <- tkframe(frameOuter,relief="groove",borderwidth=2)
	Txtscr <- tkscrollbar(tt, repeatinterval=5,
                       command=function(...)tkyview(InstrText,...))
	InstrText <- tktext(frameInstructions,bg="white",font="courier",wrap="word",yscrollcommand=function(...)tkset(Txtscr,...))
	tkgrid(InstrText,Txtscr)
	tkgrid.configure(Txtscr,sticky="ns")

	tkinsert(InstrText,"end",BFAnovaHelp()$InitInstructions)
	tkconfigure(InstrText, state="disabled")
	tkfocus(InstrText)
	
	# Bottom Frame	
	frameBottom <- tkframe(frameOuter,relief="groove",borderwidth=2)

	OnHelp=function() tkmessageBox(title="Instructions Help",message=BFAnovaHelp()$InstrHelp,icon="info",type="ok")

	OnDone=function() tclvalue(done)<-1	
	OnCancel=function() tclvalue(done)<-2
	Help.but <- tkbutton(frameBottom,text="      Help     ",command=OnHelp)
	Cancel.but <- tkbutton(frameBottom,text="      Cancel      ",command=OnCancel)
	Done.but <- tkbutton(frameBottom,text="      Next      ",command=OnDone)
	tkgrid(Help.but,column=0,row=0)
	tkgrid(Cancel.but,column=1,row=0)
	tkgrid(Done.but,column=2,row=0)
	tkgrid.configure(Done.but,sticky="w")

	tkbind(tt,"<Destroy>",OnCancel)
	tkbind(tt,"<Escape>",OnCancel)
	tkbind(tt,"<Return>",OnDone)

	#Put all frames up
	tkgrid(frameOuter)

	tkgrid(frameBottom,row=2,column=0,columnspan=2)
	tkgrid(frameInstructions,row=1,column=0,columnspan=2)
	tkwait.variable(done)

	doneVal <- as.integer(tclvalue(done))
      	tkgrab.release(tt)
	tkdestroy(tt)
	
	if(doneVal==1) retVal="OK"
	if(doneVal==2) retVal=NULL	

	return(retVal)

}



BFAnovaSelectColumns <- function(data){	
	cols=colnames(data)
	thisEnv=environment()
	SelCols=NULL
	
	tt<-tktoplevel()
	done <- tclVar(0) # Are we done yet?
	tkgrab.set(tt)
	tktitle(tt) = "Select columns of interest"
	frameOuter <- tkframe(tt)
	tkgrid(tklabel(frameOuter,text="Data Column Selection",font=BFAnovaFonts()$fontHeading),row=0)
	

	# Left Frame
	frameLeft <- tkframe(frameOuter,relief="flat",borderwidth=2)
	
	scr.AllCols <- tkscrollbar(frameLeft, repeatinterval=5,
				   command=function(...)tkyview(tl.AllCols,...))
	tl.AllCols<-tklistbox(frameLeft,height=20,selectmode="single",
				   yscrollcommand=function(...)tkset(scr.AllCols,...),background="white")
	tkgrid(tl.AllCols,scr.AllCols)
	tkgrid.configure(scr.AllCols,rowspan=20,sticky="nsw")
	tkgrid(tklabel(frameLeft,text="Available Data Columns"))
	for (i in 1:length(cols))
	{
    		tkinsert(tl.AllCols,"end",cols[i])
	}
	tkselection.set(tl.AllCols,0)

	# Middle Upper Frame
	frameMiddleUpper <- tkframe(frameOuter,relief="flat",borderwidth=2)
	
	OnDependent=function()
	{
		ColIndex <- as.integer(tkcurselection(tl.AllCols))  
		if(tclvalue(tkcurselection(tl.AllCols))!=""){
			oldVal=tclvalue(DepCol)
			ColName=cols[as.numeric(tkcurselection(tl.AllCols))+1]
			tclvalue(DepCol)<-ColName
			if(oldVal!="<not set>")
			{
				assign("cols",c(cols,oldVal),thisEnv)
				tkinsert(tl.AllCols,"end",oldVal)
			}
			assign("cols",cols[-match(ColName,cols)],thisEnv)
			tkdelete(tl.AllCols,ColIndex)
			tkselection.set(tl.AllCols,0)
		}
	}

	Dep.but <- tkbutton(frameMiddleUpper,text="Dependent Variable--->",command=OnDependent)
	tkgrid(Dep.but)


	# Middle Lower Frame
	frameMiddleLower <- tkframe(frameOuter,relief="flat",borderwidth=2)
	
	OnCateg=function()
	{
		ColIndex <- as.integer(tkcurselection(tl.AllCols))
   		if(tclvalue(tkcurselection(tl.AllCols))!=""){
   			ColName <- cols[as.numeric(tkcurselection(tl.AllCols))+1]
  	 		tkdelete(tl.AllCols,ColIndex)
			tkinsert(tl.SelCols,"end",ColName)
			assign("SelCols",rbind(SelCols,c(ColName,"CATEG")),thisEnv)
			assign("cols",cols[-match(ColName,cols)],thisEnv)
			tkselection.set(tl.AllCols,0)
		}
	}

	OnRemove=function()
	{
		ColIndex <- as.integer(tkcurselection(tl.SelCols))
   		if(tclvalue(tkcurselection(tl.SelCols))!=""){
   			if(length(SelCols)==2){
   				ColName <- SelCols[1]   		
   			}else{
   				ColName <- SelCols[as.numeric(tkcurselection(tl.SelCols))+1,1]
   			}
   			tkdelete(tl.SelCols,ColIndex)
			tkinsert(tl.AllCols,"end",ColName)
			if(length(SelCols)==2){
				assign("SelCols",NULL,thisEnv)
			}else{
				assign("SelCols",SelCols[-match(ColName,SelCols[,1]),],thisEnv)	
			}
			assign("cols",c(cols,ColName),thisEnv)	
			tkselection.set(tl.SelCols,0)
		}
	}	
	CATEG.but <- tkbutton(frameMiddleLower,text="Categorical Variable --->",command=OnCateg)
	Remove.but <- tkbutton(frameMiddleLower,text="<---Remove Variable",command=OnRemove)
	tkgrid(CATEG.but)
	tkgrid(Remove.but)
	
	# Right Upper Frame
	frameRightUpper <- tkframe(frameOuter,relief="flat",borderwidth=2)
	
	DepCol <- tclVar("<not set>")
	entry.Dep <-tkentry(frameRightUpper,width="15",textvariable=DepCol)

	tkgrid(tklabel(frameRightUpper,text="Dependent Variable:"),entry.Dep)

	tkconfigure(entry.Dep, state="readonly")


	# Right Lower Frame
	frameRightLower <- tkframe(frameOuter,relief="flat",borderwidth=2)

	
	scr.SelCols <- tkscrollbar(frameRightLower, repeatinterval=5,
				   command=function(...)tkyview(tl.SelCols,...))
	tl.SelCols<-tklistbox(frameRightLower,height=10,selectmode="single",
				   yscrollcommand=function(...)tkset(scr.SelCols,...),background="white")
	tkgrid(tl.SelCols,scr.SelCols)
	tkgrid.configure(scr.SelCols,rowspan=10,sticky="nsw")
	tkgrid(tklabel(frameRightLower,text="Selected Data Columns"))
			
	# Bottom Frame	
	frameBottom <- tkframe(frameOuter,relief="groove",borderwidth=2)

	OnHelp=function() tkmessageBox(title="Column Selection Help",message=BFAnovaHelp()$ColSelHelp,icon="info",type="ok")

	OnDone=function(){ 
		error=""
		## Sanity check on columns
		depcolnum = match(tclvalue(DepCol),colnames(data))
		if(is.na(depcolnum)){
			error=paste(error,"No dependent variable specified.\n\n",sep="")
		}
		
		if(is.null(SelCols))
			error=paste(error,"You must specify at least independent variable.\n\n",sep="")
		
		if(identical("",error)){ 
			tclvalue(done)<-1	
		}else{
			modalDialog(title="Error",message=error,reFocus=tt)
		}
	}
	OnCancel=function() tclvalue(done)<-2
	Help.but <- tkbutton(frameBottom,text="      Help     ",command=OnHelp)
	Cancel.but <- tkbutton(frameBottom,text="      Cancel      ",command=OnCancel)
	Done.but <- tkbutton(frameBottom,text="      Next      ",command=OnDone)
	tkgrid(Help.but,column=0,row=0)
	tkgrid(Cancel.but,column=1,row=0)
	tkgrid(Done.but,column=2,row=0)
	tkgrid.configure(Done.but,sticky="w")
	
	tkbind(tt,"<Destroy>",OnCancel)
	tkbind(tt,"<Escape>",OnCancel)
	tkbind(tt,"<Return>",OnDone)

	# Put all frames up
	tkgrid(frameOuter)
	HorizSep=ttkseparator(frameOuter,orient="horizontal")
	tkgrid(HorizSep,row=2,column=1,columnspan=2)
	tkgrid(frameLeft,row=1,column=0,rowspan=3)
	tkgrid(frameMiddleUpper,row=1,column=1)
	tkgrid(frameMiddleLower,row=3,column=1)
	tkgrid(frameRightUpper,row=1,column=2)
	tkgrid(frameRightLower,row=3,column=2)
	tkgrid(frameBottom,row=4,columnspan=3)
	
	tkgrid.configure(HorizSep,sticky="we")
	tkgrid.configure(frameBottom,sticky="we")
	
	tkwait.variable(done)

	doneVal <- as.integer(tclvalue(done))
      	tkgrab.release(tt)
	tkdestroy(tt)
	
	if(doneVal==1) retVal=list(dependent=tclvalue(DepCol),selectedcols=SelCols)
	if(doneVal==2) retVal=list(0)	

	return(retVal)
}

modalDialog <- function(title,message,reFocus)
{
  dlg <- tktoplevel()
  tkwm.deiconify(dlg)
  tkgrab.set(dlg)
  tkfocus(dlg)
  tkwm.title(dlg,title)
  ReturnVal <- NULL
  label1 <- tklabel(dlg,text=message)
  tkgrid(label1)

  onOK <- function()
  {
    ReturnVal <<- "OK"
    tkgrab.release(dlg)
    tkdestroy(dlg)
    tkfocus(reFocus)
  }
  OK.but     <-tkbutton(dlg,text="   OK   ",command=onOK)
  tkgrid(OK.but)

  tkfocus(dlg)
  tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(reFocus)})
   tkwait.window(dlg)

  return(ReturnVal)

}


BFAnovaSelectEffects <- function(allEffects){	
	thisEnv=environment()

	myEffs=numeric(0)


	tt<-tktoplevel()
	done <- tclVar(0) # Are we done yet?
	tkgrab.set(tt)
	tktitle(tt) = "Build Model"
	frameOuter <- tkframe(tt)
	tkgrid(tklabel(frameOuter,text="Model Building",font=BFAnovaFonts()$fontHeading),row=0,columnspan=3)

	# Left frame
	frameLeft <- tkframe(frameOuter,relief="flat",borderwidth=2)
	xScr       <- tkscrollbar(frameLeft,command=function(...)tkxview(treeWidget,...),orient="horizontal")
	yScr       <- tkscrollbar(frameLeft,command=function(...)tkyview(treeWidget,...))
	treeWidget <- tkwidget(frameLeft,"Tree",xscrollcommand=function(...)tkset(xScr,...),
                                 yscrollcommand=function(...)tkset(yScr,...),width=30,height=20)
	tkgrid(treeWidget,yScr)
	tkgrid.configure(yScr,stick="nsw")
	tkgrid(xScr)
	tkgrid.configure(xScr,stick="new")
	tkgrid(tklabel(frameLeft,text="Available"))


	for(i in 1:max(allEffects[,1])){
		if(i==1){
			nodename="Main Effects"
		}else{
			nodename=paste(i,"way")
		}
		tkinsert(treeWidget,"end","root",paste("BFEffectNode.",as.character(i),sep=""),text=nodename)
		}

	for(i in 1:dim(allEffects)[1]){
		tkinsert(treeWidget,"end",paste("BFEffectNode.",as.character(allEffects[i,1]),sep=""),as.character(allEffects[i,4]),text=as.character(allEffects[i,4]))
	}

	# Middle Upper frame
	frameMiddleUpper <- tkframe(frameOuter,relief="flat",borderwidth=2)
	
	OnAdd <- function()
	{
		choice <- tclvalue(tcl(treeWidget,"selection","get"))	
		if(substr(choice,1,1)=="{" & substr(choice,nchar(choice),nchar(choice))=="}") choice=substr(choice,2,nchar(choice)-1)
		if(is.na(pmatch("BFEffectNode",choice)) & !(choice%in%myEffs)){
			tkinsert(tl.SelEffs,"end",choice)
			myEffs=cbind(myEffs,choice)
			assign("myEffs",myEffs,pos=thisEnv)
		}
	##### add 
	}

	OnRem <- function()
	{	
		if(tclvalue(tkcurselection(tl.SelEffs))!=""){
			Index <- as.integer(tkcurselection(tl.SelEffs))
			tkdelete(tl.SelEffs,Index)
			myEffs=myEffs[-(Index+1)]
			assign("myEffs",myEffs,pos=thisEnv)
		}
		##### Rem 
	}

	
	Add.but <- tkbutton(frameMiddleUpper,text="Add to Model--->",command=OnAdd)
	Rem.but <- tkbutton(frameMiddleUpper,text="<---Remove from Model",command=OnRem)
	tkgrid(Add.but)
	tkgrid(Rem.but)


	# Right Upper frame
	frameRightUpper <- tkframe(frameOuter,relief="flat",borderwidth=2)

	scr.SelEffs <- tkscrollbar(frameRightUpper, repeatinterval=5,
				   command=function(...)tkyview(tl.SelEffs,...))
	tl.SelEffs<-tklistbox(frameRightUpper,height=5,selectmode="single",
				   yscrollcommand=function(...)tkset(scr.SelEffs,...),background="white")
	tkgrid(tl.SelEffs,scr.SelEffs)
	tkgrid.configure(scr.SelEffs,rowspan=5,sticky="nsw")
	tkgrid(tklabel(frameRightUpper,text="Model Includes"))


	
	# Bottom Frame	
	frameBottom <- tkframe(frameOuter,relief="groove",borderwidth=2)

	OnHelp=function() tkmessageBox(title="Effect Selection Help",message=BFAnovaHelp()$EffSelHelp,icon="info",type="ok")

	OnDone=function(){
		error=""
		## Sanity check on effects
		if(identical(myEffs,numeric(0)))
			error=paste(error,"You must add something to the model.\n\n",sep="")
		
		if(identical("",error)){ 
			tclvalue(done)<-1	
		}else{
			modalDialog(title="Error",message=error,reFocus=tt)
		}
	}	
	OnCancel=function() tclvalue(done)<-2
	Help.but <- tkbutton(frameBottom,text="      Help     ",command=OnHelp)
	Cancel.but <- tkbutton(frameBottom,text="      Cancel      ",command=OnCancel)
	Done.but <- tkbutton(frameBottom,text="      Next      ",command=OnDone)
	tkgrid(Help.but,column=0,row=0)
	tkgrid(Cancel.but,column=1,row=0)
	tkgrid(Done.but,column=2,row=0)
	tkgrid.configure(Done.but,sticky="w")

	tkbind(tt,"<Destroy>",OnCancel)
	tkbind(tt,"<Escape>",OnCancel)
	tkbind(tt,"<Return>",OnDone)


	#Put all frames up
	tkgrid(frameOuter)
	
	tkgrid(frameLeft,row=1,column=0,rowspan=5)
	
	tkgrid(frameMiddleUpper,row=1,column=1)
	tkgrid(frameRightUpper,row=1,column=2)
	tkgrid(frameBottom,row=6,column=0,columnspan=2)
	
	tkgrid.configure(frameBottom,sticky="we")
	
	tkwait.variable(done)

	doneVal <- as.integer(tclvalue(done))
      	tkgrab.release(tt)
	tkdestroy(tt)
	
	if(!identical(myEffs,numeric(0))) allEffects[match(myEffs,allEffects[,4]),5]=1
	
	if(doneVal==1) retVal=list(SelEffs=allEffects)
	if(doneVal==2) retVal=list(0)	

	return(retVal)
}

BFAnovaPriorSetup <- function(useA=TRUE,meanMuK=3,sdMuK=10,meanMuA=0,sdMuA=10,meanMuG=0,sdMuG=10,invGamma.a0=2,invGamma.b0=1,WishartDF=2,slopeSD=100,minWishDF=0){	
	thisEnv=environment()
	
	tt<-tktoplevel()
	done <- tclVar(0) # Are we done yet?
	tkgrab.set(tt)
	tktitle(tt) = "Specify Prior Settings"
	frameOuter <- tkframe(tt)
	tkgrid(tklabel(frameOuter,text="Prior Setup",font=BFAnovaFonts()$fontHeading),row=0,columnspan=3)


	# Prior on muK frame
	frameMuKprior <- tkframe(frameOuter,relief="groove",borderwidth=2)
	meanMuK.var <- tclVar(as.character(meanMuK))
	meanMuK.Entry <-tkentry(frameMuKprior,width="5",textvariable=meanMuK.var)
	
	sdMuK.var <- tclVar(as.character(sdMuK))
	sdMuK.Entry <-tkentry(frameMuKprior,width="5",textvariable=sdMuK.var)

	tkgrid(tklabel(frameMuKprior,text="grand mean K parameter"),row=0,column=0,columnspan=3)
	tkgrid(tklabel(frameMuKprior,text="Mean of prior on the grand mean K: "),row=1,column=0)
	tkgrid(meanMuK.Entry,row=1,column=1)
	tkgrid(tklabel(frameMuKprior,text="Standard deviation of prior on the grand mean K: "),row=2,column=0)
	tkgrid(sdMuK.Entry,row=2,column=1)
	tkgrid(tklabel(frameMuKprior,text=">0"),row=2,column=2)	


	# Prior on muA frame
	frameMuAprior <- tkframe(frameOuter,relief="groove",borderwidth=2)
	meanMuA.var <- tclVar(as.character(meanMuA))
	meanMuA.Entry <-tkentry(frameMuAprior,width="5",textvariable=meanMuA.var)
	
	sdMuA.var <- tclVar(as.character(sdMuA))
	sdMuA.Entry <-tkentry(frameMuAprior,width="5",textvariable=sdMuA.var)

	tkgrid(tklabel(frameMuAprior,text="grand mean Z parameter"),row=0,column=0,columnspan=3)
	tkgrid(tklabel(frameMuAprior,text="Mean of prior on the grand mean Z: "),row=1,column=0)
	tkgrid(meanMuA.Entry,row=1,column=1)
	tkgrid(tklabel(frameMuAprior,text="Standard deviation of prior on the grand mean Z: "),row=2,column=0)
	tkgrid(tklabel(frameMuAprior,text=">0"),row=2,column=2)	
	tkgrid(sdMuA.Entry,row=2,column=1)

	# Prior on muG frame
	frameMuGprior <- tkframe(frameOuter,relief="groove",borderwidth=2)
	meanMuG.var <- tclVar(as.character(meanMuG))
	meanMuG.Entry <-tkentry(frameMuGprior,width="5",textvariable=meanMuG.var)
	
	sdMuG.var <- tclVar(as.character(sdMuG))
	sdMuG.Entry <-tkentry(frameMuGprior,width="5",textvariable=sdMuG.var)

	tkgrid(tklabel(frameMuGprior,text="grand mean G parameter"),row=0,column=0,columnspan=3)
	tkgrid(tklabel(frameMuGprior,text="Mean of prior on the grand mean G: "),row=1,column=0)
	tkgrid(meanMuG.Entry,row=1,column=1)
	tkgrid(tklabel(frameMuGprior,text="Standard deviation of prior on the grand mean G: "),row=2,column=0)
	tkgrid(tklabel(frameMuGprior,text=">0"),row=2,column=2)	
	tkgrid(sdMuG.Entry,row=2,column=1)

	# Inverse Gamma prior setup
	frameInvGammaPrior <- tkframe(frameOuter,relief="groove",borderwidth=2)
	IGa0.var <- tclVar(as.character(invGamma.a0))
	IGa0.Entry <-tkentry(frameInvGammaPrior,width="5",textvariable=IGa0.var)
	
	IGb0.var <- tclVar(as.character(invGamma.b0))
	IGb0.Entry <-tkentry(frameInvGammaPrior,width="5",textvariable=IGb0.var)

	tkgrid(tklabel(frameInvGammaPrior,text="Inverse Gamma prior on variances"),row=0,column=0,columnspan=3)
	tkgrid(tklabel(frameInvGammaPrior,text="IG 'a' parameter: "),row=1,column=0)
	tkgrid(IGa0.Entry,row=1,column=1)
	tkgrid(tklabel(frameInvGammaPrior,text=">0"),row=1,column=2)		
	tkgrid(tklabel(frameInvGammaPrior,text="IG 'b' parameter: "),row=2,column=0)
	tkgrid(IGb0.Entry,row=2,column=1)
	tkgrid(tklabel(frameInvGammaPrior,text=">0"),row=2,column=2)		


	
	# Covariance Frame
	frameWishartPrior <- tkframe(frameOuter,relief="groove",borderwidth=2)
	WishartDF.var <- tclVar(as.character(WishartDF))
	WishartDF.Entry <-tkentry(frameWishartPrior,width="5",textvariable=WishartDF.var)
	
	tkgrid(tklabel(frameWishartPrior,text="Wishart prior on covariance Matrices"),row=0,column=0,columnspan=3)
	tkgrid(tklabel(frameWishartPrior,text="Prior scale matrix is the identity matrix."),row=1,column=0,columnspan=3)	
	tkgrid(tklabel(frameWishartPrior,text="Wishart df parameter: "),row=2,column=0)
	tkgrid(WishartDF.Entry,row=2,column=1)
	tkgrid(tklabel(frameWishartPrior,text=paste(">",minWishDF,sep="")),row=2,column=2)		
	

	# Bottom Frame	
	frameBottom <- tkframe(frameOuter,relief="groove",borderwidth=2)

	OnHelp=function() tkmessageBox(title="Prior Setup Help",message=BFAnovaHelp()$PriorSetHelp,icon="info",type="ok")

	OnDone=function(){
		error=""
		## Sanity check on prior parameters
		if(as.numeric(tclvalue(IGa0.var))<=0)
			error=paste(error,"Inverse Gamma prior must have a>0.\n\n",sep="")
		if(as.numeric(tclvalue(IGb0.var))<=0)
			error=paste(error,"Inverse Gamma prior must have b>0.\n\n",sep="")
		if(as.numeric(tclvalue(sdMuK.var))<=0)
			error=paste(error,"Standard deviation of prior on muK must have sd>0.\n\n",sep="")
		if(as.numeric(tclvalue(sdMuA.var))<=0 & useA)
			error=paste(error,"Standard deviation of prior on muA must have sd>0.\n\n",sep="")
		if(as.numeric(tclvalue(sdMuG.var))<=0)
			error=paste(error,"Standard deviation of prior on muG must have sd>0.\n\n",sep="")
		if(as.numeric(tclvalue(WishartDF.var))<=minWishDF)
			error=paste(error,"For the selected covariance matrices, df must be >",minWishDF,".\n\n",sep="")
		
		if(identical("",error)){ 
			tclvalue(done)<-1	
		}else{
			modalDialog(title="Error",message=error,reFocus=tt)
		}
	}		
	OnCancel=function() tclvalue(done)<-2
	Help.but <- tkbutton(frameBottom,text="      Help     ",command=OnHelp)
	Cancel.but <- tkbutton(frameBottom,text="      Cancel      ",command=OnCancel)
	Done.but <- tkbutton(frameBottom,text="      Next      ",command=OnDone)
	tkgrid(Help.but,column=0,row=0)
	tkgrid(Cancel.but,column=1,row=0)
	tkgrid(Done.but,column=2,row=0)
	tkgrid.configure(Done.but,sticky="w")

	tkbind(tt,"<Destroy>",OnCancel)
	tkbind(tt,"<Escape>",OnCancel)
	tkbind(tt,"<Return>",OnDone)

	#Put all frames up
	tkgrid(frameOuter)


	tkgrid(frameMuKprior,row=1,column=1)
	if(useA) tkgrid(frameMuAprior,row=2,column=1)
	tkgrid(frameMuGprior,row=3,column=1)
	tkgrid(frameBottom,row=4,column=0,columnspan=2)
	tkgrid(frameInvGammaPrior,row=1,column=0)
	#tkgrid(frameSlopePrior,row=2,column=0)
	tkgrid(frameWishartPrior,row=3,column=0)
	
	
	tkwait.variable(done)

	doneVal <- as.integer(tclvalue(done))
      	tkgrab.release(tt)
	tkdestroy(tt)
	
	if(doneVal==1) retVal=list(
				IGa0=tclvalue(IGa0.var),
				IGb0=tclvalue(IGb0.var),
				meanMuK=tclvalue(meanMuK.var),
				sdMuK=tclvalue(sdMuK.var),
				meanMuA=tclvalue(meanMuA.var),
				sdMuA=tclvalue(sdMuA.var),
				meanMuG=tclvalue(meanMuG.var),
				sdMuG=tclvalue(sdMuG.var),
				WishartDF=tclvalue(WishartDF.var)
				#slopeSD=tclvalue(slopeSD.var)
				)
	if(doneVal==2) retVal=list(0)	

	return(retVal)

}

BFAnovaMCMCSetup <- function(optimMaxit=200,nIters=1000,burnin=200,testrun=FALSE,hybridEpsilon=c(.015,.025),hybridLFsteps=80,hybridMultipoint=FALSE,hybridMultipointSize=1,hybridWeight=FALSE,hybridMPWeight=TRUE,progress=10,MetHastScale=.2,useMH="0",MetHastThin=1){	
	thisEnv=environment()
	
	tt<-tktoplevel()
	done <- tclVar(0) # Are we done yet?
	tkgrab.set(tt)
	tktitle(tt) = "Set MCMC Options"
	frameOuter <- tkframe(tt)
	tkgrid(tklabel(frameOuter,text="Set MCMC Options",font=BFAnovaFonts()$fontHeading),row=0,columnspan=3)

	# HybridMC frame
	frameHybridMC<-tkframe(frameOuter,relief="groove",borderwidth=2)
	hybridEpsilonLower.var <- tclVar(as.character(hybridEpsilon[1]))
	EpsLowEntry <-tkentry(frameHybridMC,width="6",textvariable=hybridEpsilonLower.var)
	hybridEpsilonUpper.var <- tclVar(as.character(hybridEpsilon[2]))
	EpsUppEntry <-tkentry(frameHybridMC,width="6",textvariable=hybridEpsilonUpper.var)
	hybridLF.var <- tclVar(as.character(hybridLFsteps))
	LFEntry <-tkentry(frameHybridMC,width="3",textvariable=hybridLF.var)

	OnUseMH=function(){
		if(tclvalue(useMHVal)=="0")
		{
			tkconfigure(MHscale.Entry,state="disabled")
			tkconfigure(MHthin.Entry,state="disabled")
			tkconfigure(EpsLowEntry,state="normal")
			tkconfigure(EpsUppEntry,state="normal")
			tkconfigure(LFEntry,state="normal")
			tclvalue(hybridEpsilonLower.var)=as.character(hybridEpsilon[1])
			tclvalue(hybridEpsilonUpper.var)=as.character(hybridEpsilon[2])
			tclvalue(hybridLF.var)=as.character(hybridLFsteps)
		}
		if(tclvalue(useMHVal)=="1")
		{
			tkconfigure(MHscale.Entry,state="normal")		
			tclvalue(MHscale.var)=as.character(MetHastScale)
			tkconfigure(MHthin.Entry,state="normal")		
			tclvalue(MHthin.var)=as.character(MetHastThin)		
			tkconfigure(EpsLowEntry,state="disabled")
			tkconfigure(EpsUppEntry,state="disabled")
			tkconfigure(LFEntry,state="disabled")

		}
	}
	useMHcb <- tkcheckbutton(frameHybridMC,command=OnUseMH)
	useMHVal=tclVar(useMH)
	tkconfigure(useMHcb,variable=useMHVal)		
	
	MHscale.var <- tclVar(as.character(MetHastScale))
	MHscale.Entry <-tkentry(frameHybridMC,width="5",textvariable=MHscale.var)
	
	MHthin.var <- tclVar(as.character(MetHastThin))
	MHthin.Entry <-tkentry(frameHybridMC,width="5",textvariable=MHthin.var)
	
	if(useMH=="0")
		{
			tkconfigure(MHscale.Entry,state="disabled")
			tkconfigure(EpsLowEntry,state="normal")
			tkconfigure(EpsUppEntry,state="normal")
			tkconfigure(LFEntry,state="normal")
			tclvalue(hybridEpsilonLower.var)=as.character(hybridEpsilon[1])
			tclvalue(hybridEpsilonUpper.var)=as.character(hybridEpsilon[2])
			tclvalue(hybridLF.var)=as.character(hybridLFsteps)
			tkconfigure(MHthin.Entry,state="disabled")
		
		}else{
			tkconfigure(MHscale.Entry,state="normal")
			tkconfigure(MHthin.Entry,state="normal")
			tkconfigure(EpsLowEntry,state="disabled")
			tkconfigure(EpsUppEntry,state="disabled")
			tkconfigure(LFEntry,state="disabled")

		}


	
	tkgrid(tklabel(frameHybridMC,text="Hybrid Monte Carlo Settings"),row=0,column=0,columnspan=3)
	tkgrid(tklabel(frameHybridMC,text="Lower epsilon value: "),row=1,column=0)
	tkgrid(tklabel(frameHybridMC,text="Upper epsilon value: "),row=2,column=0)
	tkgrid(EpsLowEntry,row=1,column=1)
	tkgrid(EpsUppEntry,row=2,column=1)
	tkgrid(tklabel(frameHybridMC,text="Number of Leapfrog steps: "),row=3,column=0)
	tkgrid(LFEntry,row=3,column=1)
	tkgrid(tklabel(frameHybridMC,text="Use Simple Metropolis (no Hybrid)? "),row=4,column=0)
	tkgrid(useMHcb,row=4,column=1)
	tkgrid(tklabel(frameHybridMC,text="Metropolis Candidate scale: "),row=5,column=0)	
	tkgrid(MHscale.Entry,row=5,column=1)
	tkgrid(tklabel(frameHybridMC,text="Thinning Every: "),row=6,column=0)	
	tkgrid(MHthin.Entry,row=6,column=1)


	#tkconfigure(useMPWeightcb,variable=useMPWeightVal,state="disabled")
	#tkconfigure(MPSize.Entry,state="disabled")



	# Starting Values frame	
	frameStartingVals<-tkframe(frameOuter,relief="groove",borderwidth=2)
	optimMaxIter.var <- tclVar(as.character(optimMaxit))
	optimMaxitEntry <-tkentry(frameStartingVals,width="4",textvariable=optimMaxIter.var)

	tkgrid(tklabel(frameStartingVals,text="Starting Value Options"),row=0,column=0,columnspan=3)
	tkgrid(tklabel(frameStartingVals,text="Maximum Optim Iterations: "),row=1,column=0)
	tkgrid(tklabel(frameStartingVals,text="Set to 0 to disable optim() starting values."),row=2,column=0,columnspan=3)
	tkgrid(optimMaxitEntry,row=1,column=1)

	# General Setup frame	
	frameGeneral<-tkframe(frameOuter,relief="groove",borderwidth=2)
	MCMCIterations.var <- tclVar(as.character(nIters))
	MCMCIterationsEntry <-tkentry(frameGeneral,width="7",textvariable=MCMCIterations.var)
	BurninIterations.var <- tclVar(as.character(burnin))
	BurninIterationsEntry <-tkentry(frameGeneral,width="7",textvariable=BurninIterations.var)
	Progress.var <- tclVar(as.character(progress))
	ProgressEntry <-tkentry(frameGeneral,width="7",textvariable=Progress.var)
	
	isTestcb <- tkcheckbutton(frameGeneral)
	testcb.val=tclVar("0")
	tkconfigure(isTestcb,variable=testcb.val)


	tkgrid(tklabel(frameGeneral,text="General Options"),row=0,column=0,columnspan=3)
	tkgrid(tklabel(frameGeneral,text="MCMC Iterations: "),row=1,column=0)
	tkgrid(MCMCIterationsEntry,row=1,column=1)
	tkgrid(tklabel(frameGeneral,text="Burnin Iterations: "),row=2,column=0)
	tkgrid(BurninIterationsEntry,row=2,column=1)
	tkgrid(tklabel(frameGeneral,text="Report Progress Every: "),row=3,column=0)
	tkgrid(ProgressEntry,row=3,column=1)
	tkgrid(tklabel(frameGeneral,text="Test run? "),row=4,column=0)
	tkgrid(isTestcb,row=4,column=1)
	
	
	# Bottom Frame	
	frameBottom <- tkframe(frameOuter,relief="groove",borderwidth=2)

	OnHelp=function() tkmessageBox(title="MCMC Options Help",message=BFAnovaHelp()$MCMCHelp,icon="info",type="ok")

	OnDone=function(){
		error=""
		## Sanity check on MCMC parameters
		if(as.numeric(tclvalue(hybridLF.var))<=0)
			error=paste(error,"Number of leapfrog steps must be greater than 0.\n\n",sep="")
		if(as.numeric(tclvalue(hybridEpsilonLower.var))<=0 | as.numeric(tclvalue(hybridEpsilonLower.var))>as.numeric(tclvalue(hybridEpsilonUpper.var)))
			error=paste(error,"Epsilons must both be greater than 0, and the upper value must be greater than the lower value.\n\n",sep="")
		if(as.numeric(tclvalue(optimMaxIter.var))<0)
			error=paste(error,"Number of optim iterations must be 0 or more.\n\n",sep="")
		if(as.numeric(tclvalue(BurninIterations.var))<0 | as.numeric(tclvalue(BurninIterations.var))>=as.numeric(tclvalue(MCMCIterations.var)))
			error=paste(error,"Number of burnin iterations must be 0 or more, and must be less than total MCMC iterations.\n\n",sep="")
		if(as.numeric(tclvalue(MCMCIterations.var))<=0)
			error=paste(error,"The number of MCMC iterations must be >0.\n\n",sep="")
		
		if(as.numeric(tclvalue(MHthin.var))<0 | as.numeric(tclvalue(MHthin.var))>=as.numeric(tclvalue(MCMCIterations.var)))
			error=paste(error,"Number of thinning iterations must be 0 or more, and must be less than total MCMC iterations.\n\n",sep="")
		if(as.numeric(tclvalue(MHscale.var))<=0)
			error=paste(error,"Metropolis-Hastings scale must be >0.\n\n",sep="")
		
		if(tclvalue(useMHVal)=="1" & as.numeric(tclvalue(BurninIterations.var))>=as.numeric(tclvalue(MCMCIterations.var))/as.numeric(tclvalue(MHthin.var)))
			error=paste(error,"Burnin too high - must be less than the number of samples after thinning.\n\n",sep="")
		

		if(identical("",error)){ 
			tclvalue(done)<-1	
		}else{
			modalDialog(title="Error",message=error,reFocus=tt)
		}
	}	
	OnCancel=function() tclvalue(done)<-2
	Help.but <- tkbutton(frameBottom,text="      Help     ",command=OnHelp)
	Cancel.but <- tkbutton(frameBottom,text="      Cancel      ",command=OnCancel)
	Done.but <- tkbutton(frameBottom,text="      Next      ",command=OnDone)
	tkgrid(Help.but,column=0,row=0)
	tkgrid(Cancel.but,column=1,row=0)
	tkgrid(Done.but,column=2,row=0)
	tkgrid.configure(Done.but,sticky="w")

	tkbind(tt,"<Destroy>",OnCancel)
	tkbind(tt,"<Escape>",OnCancel)
	tkbind(tt,"<Return>",OnDone)

	#Put all frames up
	tkgrid(frameOuter)

	tkgrid(frameBottom,row=2,column=0,columnspan=2)
	tkgrid(frameGeneral,row=1,column=0)
	tkgrid(frameHybridMC,row=1,column=2)	
	tkgrid(frameStartingVals,row=1,column=1)
	tkwait.variable(done)


	doneVal <- as.integer(tclvalue(done))
      	tkgrab.release(tt)
	tkdestroy(tt)
	
	if(doneVal==1) retVal=list(nIter=tclvalue(MCMCIterations.var),
					burnin=tclvalue(BurninIterations.var),
					progress=tclvalue(Progress.var),
					testRun=tclvalue(testcb.val),
					optimMaxIter=tclvalue(optimMaxIter.var),
					epsLow=tclvalue(hybridEpsilonLower.var),
					epsUpp=tclvalue(hybridEpsilonUpper.var),
					leapfrog=tclvalue(hybridLF.var),
					useMH=tclvalue(useMHVal),
					MHscale=tclvalue(MHscale.var),
					MHthin=tclvalue(MHthin.var)				
					#MPsteps=tclvalue(hybridMPsteps.var)
				)
	if(doneVal==2) retVal=list(0)	

	return(retVal)

}




BFAnovaOutputSetup <- function(){	
	thisEnv=environment()
	
	tt<-tktoplevel()
	done <- tclVar(0) # Are we done yet?
	tkgrab.set(tt)
	tktitle(tt) = "Output Setup"
	frameOuter <- tkframe(tt)
	tkgrid(tklabel(frameOuter,text="Specify Output Options",font=BFAnovaFonts()$fontHeading),row=0,columnspan=3)

	# Any Output frame	
	frameOutput<-tkframe(frameOuter,relief="groove",borderwidth=2)
	
	Outputcb <- tkcheckbutton(frameOutput)
	Outputcb.val=tclVar("1")
	tkconfigure(Outputcb,variable=Outputcb.val)
	outputlab = tklabel(frameOutput,text="Output files selected below?")
	tkgrid(outputlab,row=0,column=1)
	tkgrid(Outputcb,row=0,column=0)
	tkgrid(tklabel(frameOutput,text=paste("All files will be saved in ",getwd(),sep="")),row=1,column=1)
	tkgrid.configure(outputlab,sticky="w")

	
	# Chains frame	
	frameChains<-tkframe(frameOuter,relief="groove",borderwidth=2)
	
	ChnEffcb <- tkcheckbutton(frameChains)
	ChnEffcb.val=tclVar("1")
	tkconfigure(ChnEffcb,variable=ChnEffcb.val)
	ChnMeanscb <- tkcheckbutton(frameChains)
	ChnMeanscb.val=tclVar("1")
	tkconfigure(ChnMeanscb,variable=ChnMeanscb.val)
	ChnCovcb <- tkcheckbutton(frameChains)
	ChnCovcb.val=tclVar("1")
	tkconfigure(ChnCovcb,variable=ChnCovcb.val)
	ChnCorcb <- tkcheckbutton(frameChains)
	ChnCorcb.val=tclVar("1")
	tkconfigure(ChnCorcb,variable=ChnCorcb.val)

	tkgrid(tklabel(frameChains,text="Full Chains"),row=0,column=0)
	tkgrid(tklabel(frameChains,text="Random Effects/Slopes"),row=1,column=0)
	tkgrid(ChnEffcb,row=1,column=1)
	tkgrid(tklabel(frameChains,text="Effect/Slope Means"),row=2,column=0)
	tkgrid(ChnMeanscb,row=2,column=1)
	tkgrid(tklabel(frameChains,text="Covariance Matrices"),row=3,column=0)
	tkgrid(ChnCovcb,row=3,column=1)
	tkgrid(tklabel(frameChains,text="Correlation Matrices"),row=4,column=0)
	tkgrid(ChnCorcb,row=4,column=1)

	# Quantiles frame	
	frameQuantiles<-tkframe(frameOuter,relief="groove",borderwidth=2)
	
	QntEffcb <- tkcheckbutton(frameQuantiles)
	QntEffcb.val=tclVar("1")
	tkconfigure(QntEffcb,variable=QntEffcb.val)
	QntMeanscb <- tkcheckbutton(frameQuantiles)
	QntMeanscb.val=tclVar("1")
	tkconfigure(QntMeanscb,variable=QntMeanscb.val)
	QntCovcb <- tkcheckbutton(frameQuantiles)
	QntCovcb.val=tclVar("1")
	tkconfigure(QntCovcb,variable=QntCovcb.val)
	QntCorcb <- tkcheckbutton(frameQuantiles)
	QntCorcb.val=tclVar("1")
	tkconfigure(QntCorcb,variable=QntCorcb.val)

	tkgrid(tklabel(frameQuantiles,text="Posterior Quantiles"),row=0,column=0)
	tkgrid(tklabel(frameQuantiles,text="Random Effects/Slopes"),row=1,column=0)
	tkgrid(QntEffcb,row=1,column=1)
	tkgrid(tklabel(frameQuantiles,text="Effect/Slope Means"),row=2,column=0)
	tkgrid(ChnMeanscb,row=2,column=1)
	tkgrid(tklabel(frameQuantiles,text="Covariance Matrices"),row=3,column=0)
	tkgrid(QntCovcb,row=3,column=1)
	tkgrid(tklabel(frameQuantiles,text="Correlation Matrices"),row=4,column=0)
	tkgrid(QntCorcb,row=4,column=1)

	# Quantiles frame	
	framePostMeans<-tkframe(frameOuter,relief="groove",borderwidth=2)
	
	PmnEffcb <- tkcheckbutton(framePostMeans)
	PmnEffcb.val=tclVar("1")
	tkconfigure(PmnEffcb,variable=PmnEffcb.val)
	PmnMeanscb <- tkcheckbutton(framePostMeans)
	PmnMeanscb.val=tclVar("1")
	tkconfigure(PmnMeanscb,variable=PmnMeanscb.val)
	PmnCovcb <- tkcheckbutton(framePostMeans)
	PmnCovcb.val=tclVar("1")
	tkconfigure(PmnCovcb,variable=PmnCovcb.val)
	PmnCorcb <- tkcheckbutton(framePostMeans)
	PmnCorcb.val=tclVar("1")
	tkconfigure(PmnCorcb,variable=PmnCorcb.val)

	tkgrid(tklabel(framePostMeans,text="Posterior Means"),row=0,column=0)
	tkgrid(tklabel(framePostMeans,text="Random Effects/Slopes"),row=1,column=0)
	tkgrid(PmnEffcb,row=1,column=1)
	tkgrid(tklabel(framePostMeans,text="Effect/Slope Means"),row=2,column=0)
	tkgrid(PmnMeanscb,row=2,column=1)
	tkgrid(tklabel(framePostMeans,text="Covariance Matrices"),row=3,column=0)
	tkgrid(PmnCovcb,row=3,column=1)
	tkgrid(tklabel(framePostMeans,text="Correlation Matrices"),row=4,column=0)
	tkgrid(PmnCorcb,row=4,column=1)

	# Misc1 Frame
	frameMisc1<-tkframe(frameOuter,relief="groove",borderwidth=2)
	tkgrid(tklabel(frameMisc1,text="Miscellaneous settings"),row=0,column=1)
	
	outDatacb <- tkcheckbutton(frameMisc1)
	outDatacb.val=tclVar("1")
	tkconfigure(outDatacb,variable=outDatacb.val)
	
	tkgrid(tklabel(frameMisc1,text="Parsed Data"),row=1,column=0)
	tkgrid(outDatacb,row=1,column=1)

	outRobjcb <- tkcheckbutton(frameMisc1)
	outRobjcb.val=tclVar("1")
	tkconfigure(outRobjcb,variable=outRobjcb.val)
	
	tkgrid(tklabel(frameMisc1,text="Full R object"),row=2,column=0)
	tkgrid(outRobjcb,row=2,column=1)

	VertSep1=ttkseparator(frameMisc1,orient="vertical")
	tkgrid(VertSep1,column=2,row=1,rowspan=2)
	tkgrid.configure(VertSep1,sticky="ns")

	
	outMCMCPDFcb <- tkcheckbutton(frameMisc1)
	outMCMCPDFcb.val=tclVar("1")
	tkconfigure(outMCMCPDFcb,variable=outMCMCPDFcb.val)
	
	tkgrid(tklabel(frameMisc1,text="Chains PDF"),row=1,column=3)
	tkgrid(outMCMCPDFcb,row=1,column=4)

	outACFPDFcb <- tkcheckbutton(frameMisc1)
	outACFPDFcb.val=tclVar("1")
	tkconfigure(outACFPDFcb,variable=outACFPDFcb.val)
	
	tkgrid(tklabel(frameMisc1,text="Autocorrelation PDF"),row=2,column=3)
	tkgrid(outACFPDFcb,row=2,column=4)

	VertSep2=ttkseparator(frameMisc1,orient="vertical")
	tkgrid(VertSep2,column=5,row=1,rowspan=2)
	tkgrid.configure(VertSep2,sticky="ns")

	outCovLvlscb <- tkcheckbutton(frameMisc1)
	outCovLvlscb.val=tclVar("1")
	tkconfigure(outCovLvlscb,variable=outCovLvlscb.val)
	
	tkgrid(tklabel(frameMisc1,text="Covariance Levels"),row=1,column=6)
	tkgrid(outCovLvlscb,row=1,column=7)

	outPredcb <- tkcheckbutton(frameMisc1)
	outPredcb.val=tclVar("0")
	tkconfigure(outPredcb,variable=outPredcb.val)
	
	tkgrid(tklabel(frameMisc1,text="Predicted probabilities*"),row=2,column=6)
	tkgrid(tklabel(frameMisc1,text="*(Warning: high memory usage)"),row=3,column=6)
	tkgrid(outPredcb,row=2,column=7)

	
	# Bottom Frame	
	frameBottom <- tkframe(frameOuter,relief="groove",borderwidth=2)

	OnHelp=function() tkmessageBox(title="Output Options Help",message=BFAnovaHelp()$OutputHelp,icon="info",type="ok")

	OnDone=function() tclvalue(done)<-1	
	OnCancel=function() tclvalue(done)<-2
	Help.but <- tkbutton(frameBottom,text="      Help     ",command=OnHelp)
	Cancel.but <- tkbutton(frameBottom,text="      Cancel      ",command=OnCancel)
	Done.but <- tkbutton(frameBottom,text="      Next      ",command=OnDone)
	tkgrid(Help.but,column=0,row=0)
	tkgrid(Cancel.but,column=1,row=0)
	tkgrid(Done.but,column=2,row=0)
	tkgrid.configure(Done.but,sticky="w")

	tkbind(tt,"<Destroy>",OnCancel)
	tkbind(tt,"<Escape>",OnCancel)
	tkbind(tt,"<Return>",OnDone)

	#Put all frames up
	tkgrid(frameOuter)
	tkgrid(frameOutput,row=1,column=0,columnspan=3)
	tkgrid(frameChains,row=2,column=0,columnspan=1)
	tkgrid(frameQuantiles,row=2,column=1,columnspan=1)
	tkgrid(framePostMeans,row=2,column=2,columnspan=1)
	tkgrid(frameMisc1,row=3,column=0,columnspan=3)

	tkgrid(frameBottom,row=4,column=0,columnspan=2)

	tkwait.variable(done)

	doneVal <- as.integer(tclvalue(done))
      	tkgrab.release(tt)
	tkdestroy(tt)
	
	if(doneVal==1) retVal=list(doOutput=tclvalue(Outputcb.val),
								ChnEff=tclvalue(ChnEffcb.val),
								ChnMeans=tclvalue(ChnMeanscb.val),
								ChnCov=tclvalue(ChnCovcb.val),
								ChnCor=tclvalue(ChnCorcb.val),
								QntEff=tclvalue(QntEffcb.val),
								QntMeans=tclvalue(QntMeanscb.val),
								QntCov=tclvalue(QntCovcb.val),
								QntCor=tclvalue(QntCorcb.val),
								PmnEff=tclvalue(PmnEffcb.val),
								PmnMeans=tclvalue(PmnMeanscb.val),
								PmnCov=tclvalue(PmnCovcb.val),
								PmnCor=tclvalue(PmnCorcb.val),

								data=tclvalue(outDatacb.val),
								Robj=tclvalue(outRobjcb.val),
								MCMCPDF=tclvalue(outMCMCPDFcb.val),
								ACFPDF=tclvalue(outACFPDFcb.val),
								CovLvls=tclvalue(outCovLvlscb.val),
								Pred=tclvalue(outPredcb.val)
	)
	if(doneVal==2) retVal=list(0)	

	return(retVal)

}




BFAnovaGeneric <- function(){	
	thisEnv=environment()
	
	tt<-tktoplevel()
	done <- tclVar(0) # Are we done yet?
	tkgrab.set(tt)
	tktitle(tt) = "Generic Window"
	frameOuter <- tkframe(tt)
	tkgrid(tklabel(frameOuter,text="Generic Window Title",font=BFAnovaFonts()$fontHeading),row=0,columnspan=3)


	# Bottom Frame	
	frameBottom <- tkframe(frameOuter,relief="groove",borderwidth=2)

	OnHelp=function() tkmessageBox(title="Generic Help",message="Generic Help Content",icon="info",type="ok")

	OnDone=function() tclvalue(done)<-1	
	OnCancel=function() tclvalue(done)<-2
	Help.but <- tkbutton(frameBottom,text="      Help     ",command=OnHelp)
	Cancel.but <- tkbutton(frameBottom,text="      Cancel      ",command=OnCancel)
	Done.but <- tkbutton(frameBottom,text="      Next      ",command=OnDone)
	tkgrid(Help.but,column=0,row=0)
	tkgrid(Cancel.but,column=1,row=0)
	tkgrid(Done.but,column=2,row=0)
	tkgrid.configure(Done.but,sticky="w")

	tkbind(tt,"<Destroy>",OnCancel)
	tkbind(tt,"<Escape>",OnCancel)
	tkbind(tt,"<Return>",OnDone)

	#Put all frames up
	tkgrid(frameOuter)

	tkgrid(frameBottom,row=2,column=0,columnspan=2)

	tkwait.variable(done)

	doneVal <- as.integer(tclvalue(done))
      	tkgrab.release(tt)
	tkdestroy(tt)
	
	if(doneVal==1) retVal="OK"
	if(doneVal==2) retVal=list()	

	return(retVal)

}



