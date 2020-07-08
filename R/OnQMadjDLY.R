OnQMadjDLY<-function(){
    getfile1<-function(){
        if(!exists("ifname")) ifname<-tclvalue(tkgetOpenFile())
        else ifname<-tclvalue(tkgetOpenFile(initialfile=ifname))
        if(!nchar(ifname)){
            tkinsert(txt,"end","No file selected in QMadjDLY!\n\n")
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
    tkwm.title(tt,"QMadj_GaussianDLY")
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

    rb1 <- tkradiobutton(tt)
    rb2 <- tkradiobutton(tt)
    rb3 <- tkradiobutton(tt)
    rbValue<-tclVar("1")
    tkconfigure(rb1,variable=rbValue,value="1")
    tkconfigure(rb2,variable=rbValue,value="4")
    tkconfigure(rb3,variable=rbValue,value="12")
    tkgrid(tklabel(tt,text="Choose one of the following    "),column=1,row=4,sticky="w")
    tkgrid(tklabel(tt,text="No seasonality in trend or distribution"),column=1,row=5,sticky="w")
    tkgrid(tklabel(tt,text="Seasonal trends and distributions (DJF,MAM...)"),column=1,row=6,sticky="w")
    tkgrid(tklabel(tt,text="Monthly trends and distributions (Jan,Feb...) "),column=1,row=7,sticky="w")
    tkgrid(rb1,column=2,row=5)
    tkgrid(rb2,column=2,row=6)
    tkgrid(rb3,column=2,row=7)

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
        rbVal<-as.character(tclvalue(rbValue))
        if(oflg==0) {
            tkinsert(txt,"end",GuiErrorMSG)
            tkfocus(tt)
        }
        else{
            if(rbVal=="1"){
                itmp<-QMadj.GaussianDLY(InSeries=ifname,output=ofname,InCs=iDIpsName,
                                        MissingValueCode=MissingStr,Iadj=as.numeric(AdjStr),
                                        Mq=as.numeric(Mq0Str),Ny4a=as.numeric(Ny4aStr),GUI=TRUE)
                if(itmp<0){
                    tkinsert(txt,"end",ErrorMSG)
                    tkdestroy(tt)
                    tkfocus(main)
                    return()
                }
            }
            else{
                if(rbVal=="4"){
                    Snames<-c("DJF","MAM","JJA","SON")
                    Sstrts<-c(1201,301,601,901)
                    Sends<-c(229,531,831,1130)
                    Nseas<-4
                }
                else if(rbVal=="12"){
                    Snames<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
                    Sstrts<-c((1:12)*100+1)
                    Sends<-c(131,229,331,430,531,630,731,831,930,1031,1130,1231)
                    Nseas<-12
                }
                if(is.csv(ifname)) idata<-read.csv(ifname)
                else idata<-read.table(ifname)
                mmdd<-idata[,2]*100+idata[,3]
                for(ith in 1:Nseas){
                    if(Sstrts[ith]>Sends[ith])
                        odata<-idata[mmdd>=Sstrts[ith]|mmdd<=Sends[ith],]
                    else
                        odata<-idata[mmdd>=Sstrts[ith]&mmdd<=Sends[ith],]
                    iftmp<-paste(ofname,"_",Snames[ith],".dat",sep="")
                    oftmp<-paste(ofname,"_",Snames[ith],sep="")
                    write.table(odata,file=iftmp,col.names=F,row.names=F)
                    itmp<-QMadj.GaussianDLY(InSeries=iftmp,output=oftmp,InCs=iDIpsName,
                                            MissingValueCode=MissingStr,Iadj=as.numeric(AdjStr),
                                            Mq=as.numeric(Mq0Str),Ny4a=as.numeric(Ny4aStr),GUI=TRUE)
                    if(itmp<0){
                        tkinsert(txt,"end",ErrorMSG)
                        tkdestroy(tt)
                        tkfocus(main)
                        return()
                    }
                }
                ofileSout<-paste(ofname,"_QMadjDLYstat.txt",sep="")
                ofileAout<-paste(ofname,"_QMadjDLY.dat",sep="")
                if(file.exists(ofileSout)) file.remove(ofileSout)
                odata<-NULL
                for(ith in 1:Nseas){
                    ifileSout<-paste(ofname,"_",Snames[ith],"_QMadjDLYstat.txt",sep="")
                    ifileAout<-paste(ofname,"_",Snames[ith],"_QMadjDLY.dat",sep="")
                    iftmp<-paste(ofname,"_",Snames[ith],".dat",sep="")
                    if(ith>1) cat("\n\n",file=ofileSout,append=T)
                    cat(paste("#  ",ifname,"season:",Snames[ith],"\n"),file=ofileSout,append=T)
                    file.append(ofileSout,ifileSout)
                    file.remove(ifileSout)
                    i1data<-read.table(ifileAout)
                    odata<-rbind(odata,i1data)
                    file.remove(ifileAout)
                    file.remove(iftmp)
                }
                ymd<-odata[,1]*10000+odata[,2]*100+odata[,3]
                o1data<-odata[order(ymd),]
                write.table(o1data,file=ofileAout,col.names=F,row.names=F)
            }

            oact<-str40("ofbody",ofbody)
            ocurdir<-str40("curdir",curdir)
            ooutdir<-str40("outdir",outdir)
            if(exists("ofref")) rm("ofref",envir=.GlobalEnv)
            b20<-"                    "
            oref<-paste(b20,"                  NA",sep="")
            tkgrid(tklabel(frameMiddle,text=oact,width=40),column=2,row=5,sticky="e")
            tkgrid(tklabel(frameMiddle,text=oref,width=40),column=2,row=6,sticky="e")
            tkgrid(tklabel(frameMiddle,text=ocurdir,width=40),column=2,row=7,sticky="e")
            tkgrid(tklabel(frameMiddle,text=ooutdir,width=40),column=2,row=8,sticky="e")
            tkinsert(txt,"end","QMadjGaussianDLY finished successfully...\n")
            tkinsert(txt,"end",paste("Final output at ",outdir,"/",ofbody,"_*\n\n",sep=""))
            tkdestroy(tt)
            tkfocus(main)
            return()
        }
    }
    Ok2.but<-tkbutton(tt,text="   OK   ",command=OnOk2)
    tkgrid(Ok2.but,column=3,row=7)
    tkfocus(tt)
}

OnQMadjGaussian.wRef<-function(){
    getfile1<-function(){
        if(!exists("ifname")) ifname<-tclvalue(tkgetOpenFile())
        else ifname<-tclvalue(tkgetOpenFile(initialfile=ifname))
        if(!nchar(ifname)){
            tkinsert(txt,"end","No Base Data file selected in OnQMadjGaussian.wRef!\n\n")
            return()
        }
        otmp<-str40("ifname",ifname)
        tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=1,sticky="e")
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
        if(!exists("ifrname")) ifrname<-tclvalue(tkgetOpenFile())
        else ifrname<-tclvalue(tkgetOpenFile(initialfile=ifrname))
        if(!nchar(ifrname)){
            tkinsert(txt,"end","No Ref Data file selected in OnQMadjGaussian.wRef!\n\n")
            return()
        }
        otmp<-str40("ifrname",ifrname)
        tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=2,sticky="e")
        outdirtmp<-strsplit(ifrname,"/")[[1]]
        ofname<-outdirtmp[length(outdirtmp)]
        itmp<-strsplit(ofname,"\\.")[[1]]
        ofrbody<-itmp[1]
        assign("ifrname",ifrname,env=.GlobalEnv)
        assign("ofrbody",ofrbody,env=.GlobalEnv)
    }
    getfile3<-function(){
        if(!exists("iDIpsName")) iDIpsName<-tclvalue(tkgetOpenFile())
        else iDIpsName<-tclvalue(tkgetOpenFile(initialfile=iDIpsName))
        assign("iDIpsName",iDIpsName,env=.GlobalEnv)
        otmp<-str40("iDIpsName",iDIpsName)
        tkgrid(tklabel(tt,text=otmp,width=40),column=2,row=3,sticky="e")
    }
    tt<-tktoplevel()
    button.chg1<-tkbutton(tt,text="Change",command=getfile1)
    button.chg2<-tkbutton(tt,text="Change",command=getfile2)
    button.chg3<-tkbutton(tt,text="Change",command=getfile3)
    tkwm.title(tt,"QMadj_GaussianDLY.wRef")
    oifname<-str40("ifname",ifname)
    oifrname<-str40("ifrname",ifrname)
    oIps<-str40("iDIpsName",iDIpsName)

    fontLable<-tkfont.create(family="times",size=20,weight="bold")
    tkgrid(tklabel(tt,text="Input Base Data filename:"),column=1,row=1,sticky="w")
    tkgrid(tklabel(tt,text=oifname,width=40),column=2,row=1,sticky="e")
    tkgrid(button.chg1,column=3,row=1)

    tkgrid(tklabel(tt,text="Input Ref Data filename:"),column=1,row=2,sticky="w")
    tkgrid(tklabel(tt,text=oifrname,width=40),column=2,row=2,sticky="e")
    tkgrid(button.chg2,column=3,row=2)

    tkgrid(tklabel(tt,text="Input changepoints filename:"),column=1,row=3,sticky="w")
    tkgrid(tklabel(tt,text=oIps,width=40),column=2,row=3,sticky="e")
    tkgrid(button.chg3,column=3,row=3)

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
            #       if(!as.numeric(PlevStr)%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
            #         tkinsert(txt,"end","P_lev setting error, reset P_lev...")
            #         return()
            #       }
            ofname<-paste(paste(outdir,ofbody,sep="/"),ofrbody,sep="_")
            itmp<-QMadj.GaussianDLY.wRef(Bseries=ifname,Rseries=ifrname,InCs=iDIpsName,
                                         output=ofname,MissingValue=MissingStr,Iadj=as.numeric(AdjStr),
                                         Mq=as.numeric(Mq0Str),Ny4a=as.numeric(Ny4aStr))
            if(itmp<0){
                tkinsert(txt,"end",ErrorMSG)
                tkdestroy(tt)
                tkfocus(main)
                return()
            }
            else{
                assign("curdir",curdir,envir=.GlobalEnv)
                assign("outdir",outdir,envir=.GlobalEnv)
                assign("ifname",ifname,envir=.GlobalEnv)
                assign("ifrname",ifrname,envir=.GlobalEnv)
                assign("ofbody",ofbody,envir=.GlobalEnv)
                assign("ofname",ofname,envir=.GlobalEnv)

                oact<-str40("ofbody",ofbody)
                oref<-str40("ifrname",ifrname)
                ocurdir<-str40("curdir",curdir)
                ooutdir<-str40("outdir",outdir)
                #     UIpsName1<-paste(ofname,"_fCs.txt",sep="")
                #     UIpsName<-paste(ofname,"_mCs.txt",sep="")
                #     file.copy(from=UIpsName1,to=UIpsName,overwrite=TRUE)
                tkgrid(tklabel(frameMiddle,text=oact,width=40),column=2,row=6,sticky="e")
                tkgrid(tklabel(frameMiddle,text=oref,width=40),column=2,row=7,sticky="e")
                tkgrid(tklabel(frameMiddle,text=ocurdir,width=40),column=2,row=8,sticky="e")
                tkgrid(tklabel(frameMiddle,text=ooutdir,width=40),column=2,row=9,sticky="e")
                tkinsert(txt,"end",paste("data: ",ifname,"\noutput: ",ofname,"_*\n",sep=""))
                if(itmp<0) tkinsert(txt,"end",ErrorMSG)
                else{
                    tkinsert(txt,"end","QMadjGaussian.wRef finished successfully...\n")
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
