
Chg.Para<-function(){
  tt<-tktoplevel()
  tkwm.title(tt,"Change Parameters")

  textMissing<<-tclVar(paste(MissingStr))
  Entry.Missing<-tkentry(tt,width="10",textvariable=textMissing)
  tkgrid(tklabel(tt,text="Please enter the Missing Value Code."),sticky="w",
         column=1,row=1)
  tkgrid(Entry.Missing,column=2,row=1)

  textPlev<<-tclVar(paste(PlevStr))
  Entry.Plev<-tkentry(tt,width="10",textvariable=textPlev)
  tkgrid(tklabel(tt,text="Please enter the nominal conf. level plev value."),
         sticky="w",column=1,row=2)
  tkgrid(Entry.Plev,column=2,row=2)

  textAdj<<-tclVar(paste(AdjStr))
  Entry.Adj<-tkentry(tt,width="10",textvariable=textAdj)
  tkgrid(tklabel(tt,text="Please enter integer Iadj (0 to 10000 inclusive)"),
                 sticky="w",column=1,row=3)
  tkgrid(Entry.Adj,column=2,row=3)

  textMq0<<-tclVar(paste(Mq0Str))
  Entry.Mq0<-tkentry(tt,width="10",textvariable=textMq0)
  tkgrid(tklabel(tt,text="Please enter integer Mq (# of points for evaluating PDF)"),
         sticky="w",column=1,row=4)
  tkgrid(Entry.Mq0,column=2,row=4)

  textNy4a<<-tclVar(paste(Ny4aStr))
  Entry.Ny4a<-tkentry(tt,width="10",textvariable=textNy4a)
  tkgrid(tklabel(tt,text="Please enter integer Ny4a (>=5, or 0 for choosing the whole segment)"),
                 sticky="w",column=1,row=5)
  tkgrid(Entry.Ny4a,column=2,row=5)

  OnOk1<-function(){
    oflg<-1
    GuiErrorMSG<-NULL
    MissingStr<-tclvalue(textMissing)
    olen<-40
    assign("MissingStr",MissingStr,envir=.GlobalEnv)
    if(nchar(MissingStr)>olen) {
      GuiErrorMSG<-"MissingCode length error!\n"
      oflg<-0
    }

    PlevStr<-tclvalue(textPlev)
    if(!as.numeric(PlevStr)%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
      GuiErrorMSG<-paste(GuiErrorMSG,"plev must be one of these: 0.75,0.80,0.90,0.95,0.99,0.9999. Please re-enter\n")
      oflg<-0
    }

    AdjStr<-tclvalue(textAdj)
    if(!as.numeric(AdjStr)%in%c(0:10000)){
      GuiErrorMSG<-paste(GuiErrorMSG,"Integer Iadj must be between 0 and 10000 inclusive, please re-enter\n")
      oflg<-0
    }

    Mq0Str<-tclvalue(textMq0)
    if(!as.numeric(Mq0Str)%in%c(0:100)){
      GuiErrorMSG<-paste(GuiErrorMSG,"Mq setting must be an integer between 0 and 100, please re-enter\n")
      oflg<-0
    }

    Ny4aStr<-tclvalue(textNy4a)
    if(as.numeric(Ny4aStr)!=as.integer(Ny4aStr)){
       GuiErrorMSG<-paste(GuiErrorMSG,"Ny4a must be an integer, please re-enter\n")
       oflg<-0
    }

    if(oflg==0) {
      tkinsert(txt,"end",GuiErrorMSG)
      tkfocus(tt)
    }
    else{
      assign("MissingStr",MissingStr,envir=.GlobalEnv)
      assign("PlevStr",PlevStr,envir=.GlobalEnv)
      assign("AdjStr",AdjStr,envir=.GlobalEnv)
      assign("Mq0Str",Mq0Str,envir=.GlobalEnv)
      assign("Ny4aStr",Ny4aStr,envir=.GlobalEnv)

      tkinsert(txt,"end",paste("MissingValueCode set to:",MissingStr,"..\n",sep=" "))
      tkinsert(txt,"end",paste("The nominal level plev = ",PlevStr,"..\n",sep=" "))
      tkinsert(txt,"end",paste("Current Iadj value is",AdjStr,"..\n",sep=" "))
      tkinsert(txt,"end",paste("Current Mq value is",Mq0Str,"..\n",sep=" "))
      tkinsert(txt,"end",paste("Current Ny4a value is",Ny4aStr,"..\n",sep=" "))

      b20<-"                    "
      olen<-40
      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(MissingStr),olen)<-MissingStr
      omiss<-otmp
      tkgrid(tklabel(frameMiddle,text=omiss,width=40),column=2,row=1,sticky="e")

      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(PlevStr),olen)<-PlevStr
      oplev<-otmp
      tkgrid(tklabel(frameMiddle,text=oplev,width=40),column=2,row=2,sticky="e")

      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(AdjStr),olen)<-AdjStr
      oadj<-otmp
      tkgrid(tklabel(frameMiddle,text=oadj,width=40),column=2,row=3,sticky="e")

      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(Mq0Str),olen)<-Mq0Str
      omq0<-otmp
      tkgrid(tklabel(frameMiddle,text=omq0,width=40),column=2,row=4,sticky="e")

      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(Ny4aStr),olen)<-Ny4aStr
      ony0<-otmp
      tkgrid(tklabel(frameMiddle,text=ony0,width=40),column=2,row=5,sticky="e")

      tkdestroy(tt)
      tkfocus(main)
      return()
    }
  }

  tkbind(Entry.Ny4a,"<Return>",OnOk1)

  Ok1.but<-tkbutton(tt,text="   OK   ",command=OnOk1)
  tkbind(Entry.Missing,"<Return>",OnOk1)
  tkgrid(Ok1.but,column=1,sticky="e",row=6)
  tkfocus(tt)
}
