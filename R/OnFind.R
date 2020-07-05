
OnFindU<-function(){
  getfile1<-function(){
    if(!exists("ifname")) ifname<-tclvalue(tkgetOpenFile())
    else ifname<-tclvalue(tkgetOpenFile(initialfile=ifname))
    if(!nchar(ifname)){
      tkinsert(txt,"end","No file selected in FindUD!\n\n")
      return()
    }
    otmp<-str40("ifname",ifname)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=2,sticky="e")
    outdirtmp<-strsplit(ifname,"/")[[1]]
    if(length(outdirtmp)<=2){
      curdir<-paste(outdirtmp[1])
      outdir<-paste(outdirtmp[1],"output",sep=":/")
    }
    else{
      curdir<-outdirtmp[1]
      for(i in 2:(length(outdirtmp)-1))
        curdir<-paste(curdir,outdirtmp[i],sep="/")
        outdir<-paste(curdir,"output",sep="/")
    }
    setwd(curdir)
    if(!file.exists(outdir)) dir.create(outdir)
    ofname<-outdirtmp[length(outdirtmp)]
    itmp<-strsplit(ofname,"\\.")[[1]]
    if(length(itmp)<2) ofbody<-itmp[1]
    else{
      ofbody<-itmp[1]
      if(length(itmp)>2) for(i in 2:(length(itmp)-1))
      	ofbody<-paste(ofbody,itmp[i],sep="/")
    }
    ofname<-paste(outdir,ofbody,sep="/")
    assign("curdir",curdir,envir=.GlobalEnv)
    assign("outdir",outdir,envir=.GlobalEnv)
    assign("ifname",ifname,envir=.GlobalEnv)
    assign("ofbody",ofbody,envir=.GlobalEnv)
    assign("ofname",ofname,envir=.GlobalEnv)
  }
  tt<-tktoplevel()
  button.chg1<-tkbutton(tt,text="Change",command=getfile1)
  tkwm.title(tt,"FindU")
  oifname<-str40("ifname",ifname)
  fontLable<-tkfont.create(family="times",size=20,weight="bold")
  tkgrid(tklabel(tt,text="!!Do not",font=fontLable),sticky="e",column=1,row=1)
  tkgrid(tklabel(tt,text="choose daily precipitation data!!",
         font=fontLable),row=1,column=2)
  tkgrid(tklabel(tt,text="Input Data filename:"),column=1,row=2,sticky="w")
  tkgrid(tklabel(tt,text=oifname,width=40),column=2,row=2,sticky="e")
  tkgrid(button.chg1,column=3,row=2)

  OnOk2<-function(){
    oflg<-1
    if(!file.exists(ifname)) {
      oflg<-0
      GuiErrorMSG<-paste("Input Data file ",ifname," does not exist!\n",sep="")
    }
    if(oflg==0) {
      tkinsert(txt,"end",GuiErrorMSG)
      tkfocus(tt)
    }
    else{
      outdirtmp<-strsplit(ifname,"/")[[1]]
      if(length(outdirtmp)<=2){
        curdir<-paste(outdirtmp[1])
        outdir<-paste(outdirtmp[1],"output",sep=":/")
      }
      else{
        curdir<-outdirtmp[1]
        for(i in 2:(length(outdirtmp)-1))
          curdir<-paste(curdir,outdirtmp[i],sep="/")
          outdir<-paste(curdir,"output",sep="/")
      }
      setwd(curdir)
      if(!file.exists(outdir)) dir.create(outdir)
      ofname<-outdirtmp[length(outdirtmp)]
      itmp<-strsplit(ofname,"\\.")[[1]]
      if(length(itmp)<2) ofbody<-itmp[1]
      else{
        ofbody<-itmp[1]
        if(length(itmp)>2) for(i in 2:(length(itmp)-1))
      	  ofbody<-paste(ofbody,itmp[i],sep="/")
      }
      ofname<-paste(outdir,ofbody,sep="/")
      assign("curdir",curdir,envir=.GlobalEnv)
      assign("outdir",outdir,envir=.GlobalEnv)
      assign("ofbody",ofbody,envir=.GlobalEnv)
      assign("ofname",ofname,envir=.GlobalEnv)
      if(!as.numeric(PlevStr)%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
        tkinsert(txt,"end","P_lev setting error, reset P_lev...")
        return()
      }
      itmp<-FindU(InSeries=ifname,output=ofname,
            MissingValueCode=MissingStr,p.lev=as.numeric(PlevStr),
	    Iadj=as.numeric(AdjStr),Mq=as.numeric(Mq0Str),Ny4a=as.numeric(Ny4aStr),
	    GUI=TRUE)
      if(itmp<0){
        tkinsert(txt,"end",ErrorMSG)
        return()
      }
      else{
        UIpsName1<-paste(ofname,"_1Cs.txt",sep="")
        UIpsName<-paste(ofname,"_mCs.txt",sep="")
        file.copy(from=UIpsName1,to=UIpsName,overwrite=TRUE)
        assign("iDIpsName",UIpsName,envir=.GlobalEnv)
	oact<-str40("ofbody",ofbody)
	ocurdir<-str40("curdir",curdir)
	ooutdir<-str40("outdir",outdir)
        if(exists("ofref")) rm("ofref",envir=.GlobalEnv)
        b20<-"                    "
        oref<-paste(b20,"                  NA",sep="")
        tkgrid(tklabel(frameMiddle,text=oact,width=40),column=2,row=6,sticky="e")
        tkgrid(tklabel(frameMiddle,text=oref,width=40),column=2,row=7,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ocurdir,width=40),column=2,row=8,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ooutdir,width=40),column=2,row=9,sticky="e")
        tkinsert(txt,"end","FindU finished successfully...\n")
        tkinsert(txt,"end",paste("Modify",UIpsName,"for further calculation...\n\n"))
      }
      tkdestroy(tt)
      tkfocus(main)
      return()
    }

  }
  Ok2.but<-tkbutton(tt,text="   OK   ",command=OnOk2)
  tkgrid(Ok2.but,column=3,row=3)
  tkfocus(tt)
}

OnFindUD<-function(){
  getfile1<-function(){
    if(!exists("ifname")) ifname<-tclvalue(tkgetOpenFile())
    else ifname<-tclvalue(tkgetOpenFile(initialfile=ifname))
    if(!nchar(ifname)){
      tkinsert(txt,"end","No file selected in FindUD!\n\n")
      return()
    }
    otmp<-str40("ifname",ifname)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=2,sticky="e")
    outdirtmp<-strsplit(ifname,"/")[[1]]
    if(length(outdirtmp)<=2){
      curdir<-paste(outdirtmp[1])
      outdir<-paste(outdirtmp[1],"output",sep=":/")
    }
    else{
      curdir<-outdirtmp[1]
      for(i in 2:(length(outdirtmp)-1))
        curdir<-paste(curdir,outdirtmp[i],sep="/")
        outdir<-paste(curdir,"output",sep="/")
    }
    setwd(curdir)
    if(!file.exists(outdir)) dir.create(outdir)
    ofname<-outdirtmp[length(outdirtmp)]
    itmp<-strsplit(ofname,"\\.")[[1]]
    if(length(itmp)<2) ofbody<-itmp[1]
    else{
      ofbody<-itmp[1]
      if(length(itmp)>2) for(i in 2:(length(itmp)-1))
      	ofbody<-paste(ofbody,itmp[i],sep="/")
    }
    ofname<-paste(outdir,ofbody,sep="/")
    assign("curdir",curdir,envir=.GlobalEnv)
    assign("outdir",outdir,envir=.GlobalEnv)
    assign("ifname",ifname,envir=.GlobalEnv)
    assign("ofbody",ofbody,envir=.GlobalEnv)
    assign("ofname",ofname,envir=.GlobalEnv)
  }
  getfile2<-function(){
    if(!exists("iDIpsName")) iDIpsName<-tclvalue(tkgetOpenFile())
    else iDIpsName<-tclvalue(tkgetOpenFile(initialfile=iDIpsName))
    otmp<-str40("iDIpsName",iDIpsName)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=3,sticky="e")
    assign("iDIpsName",iDIpsName,env=.GlobalEnv)
  }
  tt<-tktoplevel()
  button.chg1<-tkbutton(tt,text="Change",command=getfile1)
  button.chg2<-tkbutton(tt,text="Change",command=getfile2)
  tkwm.title(tt,"FindUD")
  oifname<-str40("ifname",ifname)
  oIpsName<-str40("iDIpsName",iDIpsName)
  fontLable<-tkfont.create(family="times",size=20,weight="bold")
  tkgrid(tklabel(tt,text="!!Do not",font=fontLable),sticky="e",column=1,row=1)
  tkgrid(tklabel(tt,text="choose daily precipitation data!!",
         font=fontLable),row=1,column=2)
  tkgrid(tklabel(tt,text="Input Data filename:"),column=1,row=2,sticky="w")
  tkgrid(tklabel(tt,text=oifname,width=40),column=2,row=2,sticky="e")
  tkgrid(button.chg1,column=3,row=2)

  tkgrid(tklabel(tt,text="Input changepoints filename:"),column=1,row=3,sticky="w")
  tkgrid(tklabel(tt,text=oIpsName,width=40),column=2,row=3,sticky="e")
  tkgrid(button.chg2,column=3,row=3)

  OnOk2<-function(){
    oflg<-1
    if(!file.exists(ifname)) {
      oflg<-0
      GuiErrorMSG<-paste("Input Data file ",ifname," does not exist!\n",sep="")
    }
    if(!file.exists(iDIpsName)) {
      oflg<-0
      GuiErrorMSG<-paste(GuiErrorMSG,"Input changepoint file ",iDIpsName,
                   " does not exist!\n",sep="")
    }
    if(oflg==0) {
      tkinsert(txt,"end",GuiErrorMSG)
      tkfocus(tt)
    }
    else{
      if(!as.numeric(PlevStr)%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
        tkinsert(txt,"end","P_lev setting error, reset P_lev...")
        return()
      }
      outdirtmp<-strsplit(ifname,"/")[[1]]
      if(length(outdirtmp)<=2){
        curdir<-paste(outdirtmp[1])
        outdir<-paste(outdirtmp[1],"output",sep=":/")
      }
      else{
        curdir<-outdirtmp[1]
        for(i in 2:(length(outdirtmp)-1))
          curdir<-paste(curdir,outdirtmp[i],sep="/")
          outdir<-paste(curdir,"output",sep="/")
      }
      setwd(curdir)
      if(!file.exists(outdir)) dir.create(outdir)
      ofname<-outdirtmp[length(outdirtmp)]
      itmp<-strsplit(ofname,"\\.")[[1]]
      if(length(itmp)<2) ofbody<-itmp[1]
      else{
        ofbody<-itmp[1]
        if(length(itmp)>2) for(i in 2:(length(itmp)-1))
      	  ofbody<-paste(ofbody,itmp[i],sep="/")
      }
      ofname<-paste(outdir,ofbody,sep="/")
      assign("curdir",curdir,envir=.GlobalEnv)
      assign("outdir",outdir,envir=.GlobalEnv)
      assign("ofbody",ofbody,envir=.GlobalEnv)
      assign("ofname",ofname,envir=.GlobalEnv)
      itmp<-FindUD(InSeries=ifname,output=ofname,InCs=iDIpsName,
            MissingValueCode=MissingStr,p.lev=as.numeric(PlevStr),
	    Iadj=as.numeric(AdjStr),Mq=as.numeric(Mq0Str),Ny4a=as.numeric(Ny4aStr),
	    GUI=TRUE)
      if(itmp<0){
        tkinsert(txt,"end",ErrorMSG)
        return()
      }
      else{
        UDIpsName1<-paste(ofname,"_pCs.txt",sep="")
        UDIpsName<-paste(ofname,"_mCs.txt",sep="")
        file.copy(from=UDIpsName1,to=UDIpsName,overwrite=TRUE)
        assign("iDIpsName",UDIpsName,envir=.GlobalEnv)
	oact<-str40("ofbody",ofbody)
	ocurdir<-str40("curdir",curdir)
	ooutdir<-str40("outdir",outdir)
        if(exists("ofref")) rm("ofref",envir=.GlobalEnv)
        b20<-"                    "
        oref<-paste(b20,"                  NA",sep="")
        tkgrid(tklabel(frameMiddle,text=oact,width=40),column=2,row=6,sticky="e")
        tkgrid(tklabel(frameMiddle,text=oref,width=40),column=2,row=7,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ocurdir,width=40),column=2,row=8,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ooutdir,width=40),column=2,row=9,sticky="e")
        tkinsert(txt,"end","FindUD finished successfully...\n")
        tkinsert(txt,"end",paste("Modify",UDIpsName,"for further calculation...\n\n"))
      }
      tkdestroy(tt)
      tkfocus(main)
      return()
    }
  }
  Ok2.but<-tkbutton(tt,text="   OK   ",command=OnOk2)
  tkgrid(Ok2.but,column=3,row=4)
  tkfocus(tt)
}

OnStepSize<-function(){
  getfile1<-function(){
    if(!exists("ifname")) ifname<-tclvalue(tkgetOpenFile())
    else ifname<-tclvalue(tkgetOpenFile(initialfile=ifname))
    if(!nchar(ifname)){
      tkinsert(txt,"end","No file selected in StepSize!\n\n")
      return()
    }
    otmp<-str40("ifname",ifname)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=2,sticky="e")
    outdirtmp<-strsplit(ifname,"/")[[1]]
    if(length(outdirtmp)<=2){
      curdir<-paste(outdirtmp[1])
      outdir<-paste(outdirtmp[1],"output",sep=":/")
    }
    else{
      curdir<-outdirtmp[1]
      for(i in 2:(length(outdirtmp)-1))
        curdir<-paste(curdir,outdirtmp[i],sep="/")
        outdir<-paste(curdir,"output",sep="/")
    }
    setwd(curdir)
    if(!file.exists(outdir)) dir.create(outdir)
    ofname<-outdirtmp[length(outdirtmp)]
    itmp<-strsplit(ofname,"\\.")[[1]]
    if(length(itmp)<2) ofbody<-itmp[1]
    else{
      ofbody<-itmp[1]
      if(length(itmp)>2) for(i in 2:(length(itmp)-1))
      	ofbody<-paste(ofbody,itmp[i],sep="/")
    }
    ofname<-paste(outdir,ofbody,sep="/")
    assign("curdir",curdir,envir=.GlobalEnv)
    assign("outdir",outdir,envir=.GlobalEnv)
    assign("ifname",ifname,envir=.GlobalEnv)
    assign("ofbody",ofbody,envir=.GlobalEnv)
    assign("ofname",ofname,envir=.GlobalEnv)
  }
  getfile2<-function(){
    if(!exists("iDIpsName")) iDIpsName<-tclvalue(tkgetOpenFile())
    else iDIpsName<-tclvalue(tkgetOpenFile(initialfile=iDIpsName))
    otmp<-str40("iDIpsName",iDIpsName)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=3,sticky="e")
    assign("iDIpsName",iDIpsName,env=.GlobalEnv)
  }
  tt<-tktoplevel()
  button.chg1<-tkbutton(tt,text="Change",command=getfile1)
  button.chg2<-tkbutton(tt,text="Change",command=getfile2)
  tkwm.title(tt,"StepSize")
  oifname<-str40("ifname",ifname)
  oIpsName<-str40("iDIpsName",iDIpsName)
  fontLable<-tkfont.create(family="times",size=20,weight="bold")
  tkgrid(tklabel(tt,text="!!Do not",font=fontLable),sticky="e",column=1,row=1)
  tkgrid(tklabel(tt,text="choose daily precipitation data!!",
         font=fontLable),row=1,column=2)
  tkgrid(tklabel(tt,text="Input Data filename:"),column=1,row=2,sticky="w")
  tkgrid(tklabel(tt,text=oifname,width=40),column=2,row=2,sticky="e")
  tkgrid(button.chg1,column=3,row=2)

  tkgrid(tklabel(tt,text="Input changepoints filename:"),column=1,row=3,sticky="w")
  tkgrid(tklabel(tt,text=oIpsName,width=40),column=2,row=3,sticky="e")
  tkgrid(button.chg2,column=3,row=3)

  OnOk2<-function(){
    oflg<-1
    if(!file.exists(ifname)) {
      oflg<-0
      GuiErrorMSG<-paste("Input Data file ",ifname," does not exist!\n",sep="")
    }
    if(!file.exists(iDIpsName)) {
      oflg<-0
      GuiErrorMSG<-paste(GuiErrorMSG,"Input changepoint file ",iDIpsName,
                   " does not exist!\n",sep="")
    }
    if(oflg==0) {
      tkinsert(txt,"end",GuiErrorMSG)
      tkfocus(tt)
    }
    else{
      if(!as.numeric(PlevStr)%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
        tkinsert(txt,"end","P_lev setting error, reset P_lev...")
        return()
      }
      outdirtmp<-strsplit(ifname,"/")[[1]]
      if(length(outdirtmp)<=2){
        curdir<-paste(outdirtmp[1])
        outdir<-paste(outdirtmp[1],"output",sep=":/")
      }
      else{
        curdir<-outdirtmp[1]
        for(i in 2:(length(outdirtmp)-1))
          curdir<-paste(curdir,outdirtmp[i],sep="/")
          outdir<-paste(curdir,"output",sep="/")
      }
      setwd(curdir)
      if(!file.exists(outdir)) dir.create(outdir)
      ofname<-outdirtmp[length(outdirtmp)]
      itmp<-strsplit(ofname,"\\.")[[1]]
      if(length(itmp)<2) ofbody<-itmp[1]
      else{
        ofbody<-itmp[1]
        if(length(itmp)>2) for(i in 2:(length(itmp)-1))
      	  ofbody<-paste(ofbody,itmp[i],sep="/")
      }
      ofname<-paste(outdir,ofbody,sep="/")
      assign("curdir",curdir,envir=.GlobalEnv)
      assign("outdir",outdir,envir=.GlobalEnv)
      assign("ofbody",ofbody,envir=.GlobalEnv)
      assign("ofname",ofname,envir=.GlobalEnv)
      itmp<-StepSize(InSeries=ifname,output=ofname,InCs=iDIpsName,
            MissingValueCode=MissingStr,p.lev=as.numeric(PlevStr),
	    Iadj=as.numeric(AdjStr),Mq=as.numeric(Mq0Str),Ny4a=as.numeric(Ny4aStr),
	    GUI=TRUE)
      if(itmp<0){
        tkinsert(txt,"end",ErrorMSG)
        return()
      }
      else{
        UIpsName1<-paste(ofname,"_fCs.txt",sep="")
        UIpsName<-paste(ofname,"_mCs.txt",sep="")
        file.copy(from=UIpsName1,to=UIpsName,overwrite=TRUE)
        oact<-str40("ofbody",ofbody)
	ocurdir<-str40("curdir",curdir)
        ooutdir<-str40("outdir",outdir)
        if(exists("ofref")) rm("ofref",envir=.GlobalEnv)
        b20<-"                    "
        oref<-paste(b20,"                  NA",sep="")
    tkgrid(tklabel(frameMiddle,text=oact,width=40),column=2,row=6,sticky="e")
    tkgrid(tklabel(frameMiddle,text=oref,width=40),column=2,row=7,sticky="e")
    tkgrid(tklabel(frameMiddle,text=ocurdir,width=40),column=2,row=8,sticky="e")
    tkgrid(tklabel(frameMiddle,text=ooutdir,width=40),column=2,row=9,sticky="e")
        tkinsert(txt,"end","StepSize finished successfully...\n")
        tkinsert(txt,"end",paste("Final output at ",outdir,"/",ofbody,"_*\n\n",sep=""))
      }
      tkdestroy(tt)
      tkfocus(main)
      return()
    }
  }
  Ok2.but<-tkbutton(tt,text="   OK   ",command=OnOk2)
  tkgrid(Ok2.but,column=3,row=4)
  tkfocus(tt)
}

OnFindU.wRef<-function(){
  getfile1<-function(){
    if(!exists("ifname")) ifname<-tclvalue(tkgetOpenFile())
    else ifname<-tclvalue(tkgetOpenFile(initialfile=ifname))
    if(!nchar(ifname)){
      tkinsert(txt,"end","No Base file selected in FindU.wRef!\n\n")
      return()
    }
    assign("ifname",ifname,envir=.GlobalEnv)
    otmp<-str40("ifname",ifname)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=2,sticky="e")
    outdirtmp<-strsplit(ifname,"/")[[1]]
    if(length(outdirtmp)<=2){
      curdir<-paste(outdirtmp[1])
      outdir<-paste(outdirtmp[1],"output",sep=":/")
    }
    else{
      curdir<-outdirtmp[1]
      for(i in 2:(length(outdirtmp)-1))
        curdir<-paste(curdir,outdirtmp[i],sep="/")
        outdir<-paste(curdir,"output",sep="/")
    }
    setwd(curdir)
    if(!file.exists(outdir)) dir.create(outdir)
    ofname<-outdirtmp[length(outdirtmp)]
    itmp<-strsplit(ofname,"\\.")[[1]]
    if(length(itmp)<2) ofbody<-itmp[1]
    else{
      ofbody<-itmp[1]
      if(length(itmp)>2) for(i in 2:(length(itmp)-1))
      	ofbody<-paste(ofbody,itmp[i],sep="/")
    }
    assign("curdir",curdir,envir=.GlobalEnv)
    assign("outdir",outdir,envir=.GlobalEnv)
    assign("ofbody",ofbody,envir=.GlobalEnv)
  }
  getfile2<-function(){
    if(!exists("ifrname")) ifrname<-tclvalue(tkgetOpenFile())
    else ifrname<-tclvalue(tkgetOpenFile(initialfile=ifrname))
    if(!nchar(ifrname)){
      tkinsert(txt,"end","No Ref file selected in FindU.wRef!\n\n")
      return()
    }
    assign("ifrname",ifrname,env=.GlobalEnv)
    outdirtmp<-strsplit(ifrname,"/")[[1]]
    ofname<-outdirtmp[length(outdirtmp)]
    itmp<-strsplit(ofname,"\\.")[[1]]
    ofrbody<-itmp[1]
    assign("ofrbody",ofrbody,env=.GlobalEnv)
    otmp<-str40("ifrname",ifrname)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=3,sticky="e")
  }
  tt<-tktoplevel()
  button.chg1<-tkbutton(tt,text="Change",command=getfile1)
  button.chg2<-tkbutton(tt,text="Change",command=getfile2)
  tkwm.title(tt,"FindU.wRef")
  oifname<-str40("ifname",ifname)
  oifrname<-str40("ifrname",ifrname)
  
  fontLable<-tkfont.create(family="times",size=20,weight="bold")
  tkgrid(tklabel(tt,text="!!Do not",font=fontLable),sticky="e",column=1,row=1)
  tkgrid(tklabel(tt,text="choose daily precipitation data!!",
         font=fontLable),row=1,column=2)
  tkgrid(tklabel(tt,text="Input Base Data filename:"),column=1,row=2,sticky="w")
  tkgrid(tklabel(tt,text=oifname,width=40),column=2,row=2,sticky="e")
  tkgrid(button.chg1,column=3,row=2)

  tkgrid(tklabel(tt,text="Input Ref Data filename:"),column=1,row=3,sticky="w")
  tkgrid(tklabel(tt,text=oifrname,width=40),column=2,row=3,sticky="e")
  tkgrid(button.chg2,column=3,row=3)

  OnOk3<-function(){
    oflg<-1
    if(!file.exists(ifname)) {
      oflg<-0
      GuiErrorMSG<-paste("Input Data file ",ifname," does not exist!\n",sep="")
    }
    if(!file.exists(ifrname)) {
      oflg<-0
      GuiErrorMSG<-paste("Input Ref file ",ifrname," does not exist!\n",sep="")
    }
    if(oflg==0) {
      tkinsert(txt,"end",GuiErrorMSG)
      tkfocus(tt)
    }
    else{
      if(!as.numeric(PlevStr)%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
        tkinsert(txt,"end","P_lev setting error, reset P_lev...")
	return()
      }
      ofname<-paste(paste(outdir,ofbody,sep="/"),ofrbody,sep="_")
      itmp<-FindU.wRef(Bseries=ifname,Rseries=ifrname,output=ofname,
            MissingValueCode=MissingStr,p.lev=as.numeric(PlevStr),
    	    Iadj=as.numeric(AdjStr),Mq=as.numeric(Mq0Str),Ny4a=as.numeric(Ny4aStr),
	    GUI=T)
      if(itmp<0){ # Error happens
        tkinsert(txt,"end",ErrorMSG)
        return()
      }
      else{
        assign("curdir",curdir,envir=.GlobalEnv)
        assign("outdir",outdir,envir=.GlobalEnv)
        assign("ifname",ifname,envir=.GlobalEnv)
        assign("ifrname",ifrname,envir=.GlobalEnv)
        assign("ofbody",ofbody,envir=.GlobalEnv)
        assign("ofname",ofname,envir=.GlobalEnv)
        UIpsName1<-paste(ofname,"_1Cs.txt",sep="")
        UIpsName<-paste(ofname,"_mCs.txt",sep="")
        file.copy(from=UIpsName1,to=UIpsName,overwrite=TRUE)
        assign("iIpsName",UIpsName1,envir=.GlobalEnv)
        assign("iDIpsName",UIpsName,envir=.GlobalEnv)
        oact<-str40("ofbody",ofbody)
        oref<-str40("ifrname",ifrname)
        ocurdir<-str40("curdir",curdir)
        ooutdir<-str40("outdir",outdir)
        tkgrid(tklabel(frameMiddle,text=oact,width=40),column=2,row=6,sticky="e")
        tkgrid(tklabel(frameMiddle,text=oref,width=40),column=2,row=7,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ocurdir,width=40),column=2,row=8,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ooutdir,width=40),column=2,row=9,sticky="e")
        tkinsert(txt,"end",paste("data: ",ifname,"\noutput: ",ofname,"_*\n",sep=""))
        if(itmp<0) tkinsert(txt,"end",ErrorMSG)
        else{
          tkinsert(txt,"end","FindU.wRef finished successfully...\n")
          tkinsert(txt,"end",paste("Modify",UIpsName,"for further calculation...\n\n"))
        }
      }
      tkdestroy(tt)
      tkfocus(main)
      return()
    }
  }
  Ok3.but<-tkbutton(tt,text="   OK   ",command=OnOk3)
  tkgrid(Ok3.but,column=3,row=4)
  tkfocus(tt)
}

OnFindUD.wRef<-function(){
  getfile1<-function(){
    if(!exists("ifname")) ifname<-tclvalue(tkgetOpenFile())
    else ifname<-tclvalue(tkgetOpenFile(initialfile=ifname))
    if(!nchar(ifname)){
      tkinsert(txt,"end","No Base Data file selected in FindU.wRef!\n\n")
      return()
    }
    assign("ifname",ifname,envir=.GlobalEnv)
    otmp<-str40("ifname",ifname)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=2,sticky="e")
    outdirtmp<-strsplit(ifname,"/")[[1]]
    if(length(outdirtmp)<=2){
      curdir<-paste(outdirtmp[1])
      outdir<-paste(outdirtmp[1],"output",sep=":/")
    }
    else{
      curdir<-outdirtmp[1]
      for(i in 2:(length(outdirtmp)-1))
        curdir<-paste(curdir,outdirtmp[i],sep="/")
        outdir<-paste(curdir,"output",sep="/")
    }
    setwd(curdir)
    if(!file.exists(outdir)) dir.create(outdir)
    ofname<-outdirtmp[length(outdirtmp)]
    itmp<-strsplit(ofname,"\\.")[[1]]
    if(length(itmp)<2) ofbody<-itmp[1]
    else{
      ofbody<-itmp[1]
      if(length(itmp)>2) for(i in 2:(length(itmp)-1))
      	ofbody<-paste(ofbody,itmp[i],sep="/")
    }
    ofname<-paste(outdir,ofbody,sep="/")
    assign("curdir",curdir,envir=.GlobalEnv)
    assign("outdir",outdir,envir=.GlobalEnv)
    assign("ofbody",ofbody,envir=.GlobalEnv)
  }
  getfile2<-function(){
    if(!exists("ifrname")) ifrname<-tclvalue(tkgetOpenFile())
    else ifrname<-tclvalue(tkgetOpenFile(initialfile=ifrname))
    if(!nchar(ifrname)){
      tkinsert(txt,"end","No Ref Data file selected in FindU.wRef!\n\n")
      return()
    }
    assign("ifrname",ifrname,env=.GlobalEnv)
    outdirtmp<-strsplit(ifrname,"/")[[1]]
    ofname<-outdirtmp[length(outdirtmp)]
    itmp<-strsplit(ofname,"\\.")[[1]]
    ofrbody<-itmp[1]
    assign("ofrbody",ofrbody,env=.GlobalEnv)
    otmp<-str40("ifrname",ifrname)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=3,sticky="e")
  }
  getfile3<-function(){
    if(!exists("iDIpsName")) iDIpsName<-tclvalue(tkgetOpenFile())
    else iDIpsName<-tclvalue(tkgetOpenFile(initialfile=iDIpsName))
    assign("iDIpsName",iDIpsName,env=.GlobalEnv)
    otmp<-str40("iDIpsName",iDIpsName)
    tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=4,sticky="e")
  }
  tt<-tktoplevel()
  button.chg1<-tkbutton(tt,text="Change",command=getfile1)
  button.chg2<-tkbutton(tt,text="Change",command=getfile2)
  button.chg3<-tkbutton(tt,text="Change",command=getfile3)
  tkwm.title(tt,"FindUD.wRef")
  oifname<-str40("ifname",ifname)
  oifrname<-str40("ifrname",ifrname)
  oIps<-str40("iDIpsName",iDIpsName)
  
  fontLable<-tkfont.create(family="times",size=20,weight="bold")
  tkgrid(tklabel(tt,text="!!Do not",font=fontLable),sticky="e",column=1,row=1)
  tkgrid(tklabel(tt,text="choose daily precipitation data!!",
         font=fontLable),row=1,column=2)
  tkgrid(tklabel(tt,text="Input Base Data filename:"),column=1,row=2,sticky="w")
  tkgrid(tklabel(tt,text=oifname,width=40),column=2,row=2,sticky="e")
  tkgrid(button.chg1,column=3,row=2)

  tkgrid(tklabel(tt,text="Input Ref Data filename:"),column=1,row=3,sticky="w")
  tkgrid(tklabel(tt,text=oifrname,width=40),column=2,row=3,sticky="e")
  tkgrid(button.chg2,column=3,row=3)

  tkgrid(tklabel(tt,text="Input changepoints filename:"),column=1,row=4,sticky="w")
  tkgrid(tklabel(tt,text=oIps,width=40),column=2,row=4,sticky="e")
  tkgrid(button.chg3,column=3,row=4)

  OnOk3<-function(){
    oflg<-1
    if(!file.exists(ifname)) {
      oflg<-0
      GuiErrorMSG<-paste("Input Base Data file ",ifname," does not exist!\n",sep="")
    }
    if(!file.exists(ifrname)) {
      oflg<-0
      GuiErrorMSG<-paste("Input Ref Data file ",ifrname," does not exist!\n",sep="")
    }
    if(!file.exists(iDIpsName)) {
      oflg<-0
      GuiErrorMSG<-paste("Input changepoint file ",iDIpsName," does not exist!\n",sep="")
    }
    if(oflg==0) {
      tkinsert(txt,"end",GuiErrorMSG)
      tkfocus(tt)
    }
    else{
      if(!as.numeric(PlevStr)%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
        tkinsert(txt,"end","P_lev setting error, reset P_lev...")
	return()
      }
      ofname<-paste(paste(outdir,ofbody,sep="/"),ofrbody,sep="_")
      itmp<-FindUD.wRef(Bseries=ifname,Rseries=ifrname,InCs=iDIpsName,
            output=ofname,MissingValueCode=MissingStr,p.lev=as.numeric(PlevStr),
    	    Iadj=as.numeric(AdjStr),Mq=as.numeric(Mq0Str),Ny4a=as.numeric(Ny4aStr),GUI=T)
      if(itmp<0){ # Error happens
        tkinsert(txt,"end",ErrorMSG)
        return()
      }
      else{
        assign("curdir",curdir,envir=.GlobalEnv)
        assign("outdir",outdir,envir=.GlobalEnv)
        assign("ifname",ifname,envir=.GlobalEnv)
        assign("ifrname",ifrname,envir=.GlobalEnv)
        assign("ofbody",ofbody,envir=.GlobalEnv)
        assign("ofname",ofname,envir=.GlobalEnv)
        UIpsName1<-paste(ofname,"_pCs.txt",sep="")
        UIpsName<-paste(ofname,"_mCs.txt",sep="")
        file.copy(from=UIpsName1,to=UIpsName,overwrite=TRUE)
        assign("iIpsName",UIpsName,envir=.GlobalEnv)
        oact<-str40("ofbody",ofbody)
        oref<-str40("ifrname",ifrname)
        ocurdir<-str40("curdir",curdir)
        ooutdir<-str40("outdir",outdir)
        tkgrid(tklabel(frameMiddle,text=oact,width=40),column=2,row=6,sticky="e")
        tkgrid(tklabel(frameMiddle,text=oref,width=40),column=2,row=7,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ocurdir,width=40),column=2,row=8,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ooutdir,width=40),column=2,row=9,sticky="e")
        tkinsert(txt,"end",paste("data: ",ifname,"\noutput: ",ofname,"_*\n",sep=""))
        if(itmp<0) tkinsert(txt,"end",ErrorMSG)
        else{
          tkinsert(txt,"end","FindUD.wRef finished successfully...\n")
          tkinsert(txt,"end",paste("Modify",UIpsName,"for further calculation...\n\n"))
        }
      }
      tkdestroy(tt)
      tkfocus(main)
      return()
    }
  }
  Ok3.but<-tkbutton(tt,text="   OK   ",command=OnOk3)
  tkgrid(Ok3.but,column=3,row=5)
  tkfocus(tt)
}
