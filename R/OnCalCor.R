
OnCalCor<-function(){
  CalCor<-function(bfname,rfname,MissingValueCode){
    ocor<-ocnt<-NA; ocomm<-''
    if(!file.exists(bfname)) {
       ocomm<<-paste("Input basefile",bfname,"does not exist!")
       return(c(ocor,ocnt,ocomm))
    }
    if(!file.exists(rfname)) {
       ocomm<-paste("Input ref file",rfname,"does not exist!")
       return(c(ocor,ocnt,ocomm))
    }
    if(is.csv(bfname)){
      itmp<-try(read.table(bfname,sep=",",header=F,na.strings=MissingValueCode,
              colClasses="real"),silent=T)
      if(inherits(itmp,"try-error")){
        ocomm<-geterrmessage()
        return(c(ocor,ocnt,ocomm))
      }
      else itable<-itmp
    }
    else{
      itmp<-try(read.table(bfname,sep="",header=F,na.strings=MissingValueCode,
            colClasses="real"),silent=T)
      if(inherits(itmp,"try-error")){
        ocomm<-geterrmessage()
        return(c(ocor,ocnt,ocomm))
      }
      else itable<-itmp
    }

    if(is.csv(rfname)){
      itmp<-try(read.table(rfname,sep=",",header=F,na.strings=MissingValueCode,
              colClasses="real"),silent=T)
      if(inherits(itmp,"try-error")){
        ocomm<-geterrmessage()
        return(c(ocor,ocnt,ocomm))
      }
      else rtable<-itmp
    }
    else{
      itmp<-try(read.table(rfname,sep="",header=F,na.strings=MissingValueCode,
            colClasses="real"),silent=T)
      if(inherits(itmp,"try-error")){
        ocomm<-geterrmessage()
        return(c(ocor,ocnt,ocomm))
      }
      else rtable<-itmp
    }
    Icy<-sort(unique(itable[,2]*100+itable[,3]))
    Nt<-length(Icy)
    yrange=range(c(itable[,1],rtable[,1]))
    oymd=NULL
    if(Nt==1){ # annual data
      oys=seq(yrange[1],yrange[2])
      odata=matrix(NA,length(oys),2)
      rownames(odata)=as.character(oys)
      odata[as.character(itable[,1]),1]=itable[,4]
      odata[as.character(rtable[,1]),2]=rtable[,4]
      pyr=oys
    }
    else if(Nt<=12){ # monthly/seasonal data
      oym=sort(rep(yrange[1]:yrange[2],12))*100+rep(1:12,length(yrange[1]:yrange[2]))
      ym.b=itable[,1]*100+itable[,2]
      ym.r=rtable[,1]*100+rtable[,2]
      odata=matrix(NA,length(oym),2)
      rownames(odata)=as.character(oym)
      odata[as.character(ym.b),1]=itable[,4]
      odata[as.character(ym.r),2]=rtable[,4]
      pyr=floor(oym/100)
    }
    else{ # count as daily
      mdays=c(31,28,31,30,31,30,31,31,30,31,30,31)
      oymd=NULL
      for(yr in yrange[1]:yrange[2]) for(mon in 1:12){
        if(mon==2&yr%%4==0&(yr%%100!=0|yr%%400==0)) dd=29 else dd=mdays[mon]
        oymd=c(oymd,yr*10000+mon*100+1:dd)
      }
      ymd.b=itable[,1]*10000+itable[,2]*100+itable[,3]
      ymd.r=rtable[,1]*10000+rtable[,2]*100+rtable[,3]
      odata=matrix(NA,length(oymd),2)
      rownames(odata)=as.character(oymd)
      odata[as.character(ymd.b),1]=itable[,4]
      odata[as.character(ymd.r),2]=rtable[,4]
      pyr=floor(oymd/10000)
    }

    data.diff=na.omit(cbind(odata[-1,1]-odata[-nrow(odata),1],odata[-1,2]-odata[-nrow(odata),2]))
    if(nrow(data.diff)>3){
      ocnt=nrow(data.diff)
      ocor=cor(data.diff[,1],data.diff[,2])*100

      labs=if(diff(yrange)>10) unique(pyr/10) else unique(pyr)
      ats=match(labs,pyr)
      ylow=min(odata,na.rm=T); yup=max(odata,na.rm=T); ylim=c(ylow,yup+(yup-ylow)*.2)
      plot(1:nrow(odata),odata[,1],type='l',col='red',ylim=ylim,xaxt='n',cex=.5,xlab=paste('Cor=',round(ocor,2)),ylab='')
#          main=basename(rfname),col.main='blue')
      lines(1:nrow(odata),odata[,2],type='l',col='blue',ylim=ylim,xaxt='n',cex=.5,xlab=paste('Cor=',round(ocor,2)),ylab='')
      title(bquote(.(paste(basename(bfname),'     '))*phantom(.(basename(rfname)))),col.main='red')
      title(bquote(phantom(.(paste(basename(bfname),'     ')))*.(basename(rfname))),col.main='blue')
    }
    else{
      ocnt=nrow(data.diff)
      ocomm='no sufficient data'
    }

    return(list(ocnt=ocnt,ocor=ocor,ocomm=ocomm))
  }

  getfile1<-function(){
    if(!exists("ifname")) ifname<-tclvalue(tkgetOpenFile())
    else ifname<-tclvalue(tkgetOpenFile(initialfile=ifname))
    if(!nchar(ifname)){
      tkinsert(txt,"end","No Base Data file selected in CalCor!\n\n")
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
    name <- tk_choose.files(caption = "Open", filters = matrix(c("All Readable Files", "*", "All Files", "*"), ncol = 2, byrow = TRUE))
    if(length(name) < 1) {
      tkinsert(txt,"end","No Base Data file selected in CalCor!\n\n")
      return()
    } else {
    ## Multiple file selected, store list, and open first dataset only
      assign("filelist", name, envir = .GlobalEnv)
      otmp<-str40("filelist",name[1])
      tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=3,sticky="e")
    }
  }

  tt<-tktoplevel()
  button.chg1<-tkbutton(tt,text="Change",command=getfile1)
  button.chg2<-tkbutton(tt,text="Change",command=getfile2)
  tkwm.title(tt,"CalCor")
  oifname<-str40("ifname",ifname)
  oifrname<-str40("ifrname",ifrname)

  fontLable<-tkfont.create(family="times",size=20,weight="bold")
  tkgrid(tklabel(tt,text="Calculating Correlations",font=fontLable),sticky="e",column=1,row=1)
# tkgrid(tklabel(tt,text="choose daily precipitation data!!",
#        font=fontLable),row=1,column=2)
  tkgrid(tklabel(tt,text="Input Base Data filename:"),column=1,row=2,sticky="w")
  tkgrid(tklabel(tt,text=oifname,width=40),column=2,row=2,sticky="e")
  tkgrid(button.chg1,column=3,row=2)

  tkgrid(tklabel(tt,text="Input Ref Data filename:(multiple files)"),column=1,row=3,sticky="w")
  tkgrid(tklabel(tt,text=oifrname,width=40),column=2,row=3,sticky="e")
  tkgrid(button.chg2,column=3,row=3)

  OnOk3<-function(){
    oflg<-1
    if(!file.exists(ifname)) {
      oflg<-0
      GuiErrorMSG<-paste("Input Base Data file ",ifname," does not exist!\n",sep="")
    }
    if(oflg==0) {
      tkinsert(txt,"end",GuiErrorMSG)
      tkfocus(tt)
    }

    ofilePdf=paste(outdir,'/',ofbody,'_Cor.pdf',sep='')
    oflog=paste(outdir,'/',ofbody,'_Cor.log',sep='')

    pdf(file=ofilePdf,onefile=T,paper='letter')
    op <- par(no.readonly = TRUE) # the whole list of settable par's.
    par(mfrow=c(2,1))
    par(mar=c(3,4,3,2)+.1,cex.main=.8,cex.lab=.8,cex.axis=.8,cex=.8)

    for(rfname in filelist) {
      otmp=CalCor(ifname,rfname,MissingStr)
      cat(paste(basename(ifname),basename(rfname),'Cor=',round(otmp$ocor,2),
                ', total sample number:',otmp$ocnt,', comment:',otmp$ocomm,'\n'),
          file=oflog,append=!(rfname==filelist[1]))
    }

    par(op)
    dev.off()

    tkinsert(txt,"end","CalCor finished successfully...\n")
    tkdestroy(tt)
    tkfocus(main)
    return()
  }

  Ok3.but<-tkbutton(tt,text="   OK   ",command=OnOk3)
  tkgrid(Ok3.but,column=3,row=5)
  tkfocus(tt)
}
