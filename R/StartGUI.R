#' StartGUI
#'
#' @export
#' @import tcltk
StartGUI <- function(){
  UA1 <- 'RClimDex and RHtests software packages (all versions included), herein after called "The Library" \n'
  UA2 <- '©  Copyright, Environment Canada, 2012, \n'
  UA3 <- "The Library was created by the Climate Research Division of Environment Canada and as such all intellectual property rights (including copyright) that may exist in The Library are owned by the Government of Canada.  The Library software code is provided free under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, version 3.0 of the License. It is distributed under the terms of this license 'as-is' and has not been designed or prepared to meet any Licensee's particular requirements. Environment Canada makes no warranty, either express or implied, including but not limited to, warranties of merchantability or fitness for a particular purpose. In no event will Environment Canada be liable for any indirect, special, consequential or other damages attributed to the Licensee's use of The Library. In downloading The Library you understand and agree to these terms and those of the associated LGP License. See the GNU Lesser General Public License for more details.\n"
  UA4 <- "You should have received a copy of the GNU Lesser General Public License along with This Library; if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA."
  UserAgrement  <-  paste(UA1,UA2,UA3,UA4,sep='\n')
  
  yesAgree <- function(){
    tkdestroy(topFrame)
    GUI()
  }

  noAgree <- function(){
    tkdestroy(topFrame)
  }

  topFrame  <-  tktoplevel(pady=2, bg="gray94")
  outFrame  <-  tkframe(topFrame, padx=2, bg="gray94")
  buttonFrame  <-  tkframe(topFrame, padx=2, bg="gray94")


  tkwm.title(topFrame, "RHtestsV4 User Agreement")
  tkwm.resizable(topFrame, 0, 0)
  tkgrid(outFrame)
  tkgrid(buttonFrame)
  tkgrid.configure(outFrame, sticky="nsew")
  scrollBar  <-  tkscrollbar(outFrame, repeatinterval=5, command=function(...)tkyview(textFrame,...))
  textFrame  <-  tktext(outFrame,bg="gray94",font="courier",yscrollcommand=function(...)tkset(scrollBar,...))
  tkgrid(textFrame, scrollBar)
  tkgrid.configure(scrollBar,sticky="ns")

  tkinsert(textFrame,"end",paste(UserAgrement,sep=""))
  tkconfigure(textFrame, state="disabled")
  
  yesButton  <-  tkbutton(buttonFrame, text="I Agree", command=function()yesAgree())
  noButton  <-  tkbutton(buttonFrame, text="I Do Not Agree", command=function()noAgree())
  tkgrid(yesButton,noButton)

  tkbind(topFrame,"<Destroy>",function()noAgree())
  tkfocus(textFrame)

  # tkdestroy(topFrame)
  # GUI()
}

GUI <- function(){
  if(!exists("MissingStr")) MissingStr <- "-99.9"
  if(!exists("PlevStr")) PlevStr <- "0.95"
  if(!exists("AdjStr")) AdjStr <- "10000"
  if(!exists("Mq0Str")) Mq0Str <- "12"
  if(!exists("Ny4aStr")) Ny4aStr <- "0"
  assign("MissingStr",MissingStr,envir=.GlobalEnv)
  assign("PlevStr",PlevStr,envir=.GlobalEnv)
  assign("AdjStr",AdjStr,envir=.GlobalEnv)
  assign("Mq0Str",Mq0Str,envir=.GlobalEnv)
  assign("Ny4aStr",Ny4aStr,envir=.GlobalEnv)
  main <- tktoplevel()
  assign("main",main,envir=.GlobalEnv)
  tkwm.title(main,"RHtestsV4")
  fontHeading <- tkfont.create(family="times",size=40,weight="bold", slant="italic")
  fontLable <- tkfont.create(family="times",size=15,weight="bold")

  frameOverall <- tkframe(main)
  frameUpper <- tkframe(frameOverall,relief="groove",borderwidth=2)
  assign("frameUpper",frameUpper,env=.GlobalEnv)
  frameMiddle <- tkframe(frameOverall,relief="groove",borderwidth=2)
  assign("frameMiddle",frameMiddle,env=.GlobalEnv)
  tkgrid(tklabel(frameUpper,text="RHtests V4",font=fontHeading),column=2,row=1)

  ChgPara.but <- tkbutton(frameUpper,text="Change Pars",width=14,command=Chg.Para)
  Rfile.but <- tkbutton(frameUpper,text="Transform Data",width=14,command=Read.file)
  FindU.but <- tkbutton(frameUpper,text= "FindU",width=14,command=OnFindU)
  
  FindUD.but <- tkbutton(frameUpper,text="FindUD",width=14,command=OnFindUD)
  StepSize.but <- tkbutton(frameUpper,text="StepSize",width=14,command=OnStepSize)
  Cancel.but <- tkbutton(frameUpper,text="Quit",width=14,command=OnQuit)
  FindU.wRef.but <- tkbutton(frameUpper,text= "FindU.wRef",width=14,
                            command=OnFindU.wRef)
  FindUD.wRef.but <- tkbutton(frameUpper,text="FindUD.wRef",width=14,
                            command=OnFindUD.wRef)
  StepSize.wRef.but <- tkbutton(frameUpper,text="StepSize.wRef",width=14,
                              command=OnStepSize.wRef)
  CalCor.but <- tkbutton(frameUpper,text='CalCol',width=14,command=OnCalCor)
  QMadjDLY.but <- tkbutton(frameUpper,text="QMadj",width=14,
                          command=OnQMadjDLY)
  QMadjDLY.wRef.but <- tkbutton(frameUpper,text="QMadj.wRef",width=14,
                              command=OnQMadjGaussian.wRef)
  # FindCor.but <- tkbutton(frameUpper,text='FindCor',width=14,command=OnCalCor)
  olen <- 40
  b20 <- "                    "
  if(nchar(MissingStr)>olen) stop("MissingValueCode length error!")
  otmp <- paste(b20,b20,sep="")
  substr(otmp,olen+1-nchar(MissingStr),olen) <- MissingStr
  omiss <- otmp
  if(nchar(PlevStr)>olen) stop("P_lev length error!")
  otmp <- paste(b20,b20,sep="")
  substr(otmp,olen+1-nchar(PlevStr),olen) <- PlevStr
  oplev <- otmp
  otmp <- paste(b20,b20,sep="")
  substr(otmp,olen+1-nchar(AdjStr),olen) <- AdjStr
  oadj <- otmp
  otmp <- paste(b20,b20,sep="")
  substr(otmp,olen+1-nchar(Mq0Str),olen) <- Mq0Str
  omq0 <- otmp
  otmp <- paste(b20,b20,sep="")
  substr(otmp,olen+1-nchar(Ny4aStr),olen) <- Ny4aStr
  oNy4a <- otmp
  if(!exists("ofbody")) otmp <- paste(b20,"                  NA",sep="")
  else{
    if(nchar(ofbody)>olen) otmp <- paste("...",substr(ofbody,nchar(ofbody)-olen+4,nchar(ofbody)),sep="")
    else{
      otmp <- paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(ofbody),olen) <- ofbody
    }
  }
  oact <- otmp
  if(!exists("ofref")) otmp <- paste(b20,"                  NA",sep="")
  else{
    if(nchar(ofref)>olen) otmp <- paste("...",substr(ofref,nchar(ofref)-olen+4,nchar(ofref)),sep="")
    else{
      otmp <- paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(ofref),olen) <- ofref
    }
  }
  oref <- otmp
  if(!exists("curdir")) otmp <- paste(b20,"                  NA",sep="")
  else{
    if(nchar(curdir)>olen) otmp <- paste("...",substr(curdir,nchar(curdir)-olen+4,nchar(curdir)),sep="")
    else{
      otmp <- paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(curdir),olen) <- curdir
    }
  }
  ocurdir <- otmp
  if(!exists("outdir")) otmp <- paste(b20,"                  NA",sep="")
  else{
    if(nchar(outdir)>olen) otmp <- paste("...",substr(outdir,nchar(outdir)-olen+4,nchar(outdir)),sep="")
    else{
      otmp <- paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(outdir),olen) <- outdir
    }
  }
  ooutdir <- otmp
  # arrange menu buttons and labels
  tkgrid(tklabel(frameUpper,text="          "),Rfile.but,ChgPara.but,Cancel.but)
  tkgrid(tklabel(frameUpper,text="PMT and t tests:",font=fontLable),sticky="w")
  tkgrid(FindU.wRef.but,column=1,row=3)
  tkgrid(FindUD.wRef.but,column=2,row=3)
  tkgrid(StepSize.wRef.but,column=3,row=3)
  tkgrid(tklabel(frameUpper,text="PMF and F tests:",font=fontLable),sticky="w")
  tkgrid(FindU.but,column=1,row=4)
  tkgrid(FindUD.but,column=2,row=4)
  tkgrid(StepSize.but,column=3,row=4)
  tkgrid(tklabel(frameUpper,text="  To adjust daily",font=fontLable),sticky="w")
  tkgrid(tklabel(frameUpper,text="Gaussian data:",font=fontLable),sticky="w",column=1,row=5)
  tkgrid(tklabel(frameUpper,text="Calculate Correlation", font=fontLable),CalCor.but)
  # tkgrid(CalCor.but,column=2,row=6)
  tkgrid(QMadjDLY.wRef.but,column=2,row=5)
  tkgrid(QMadjDLY.but,column=3,row=5)
  tkgrid(tklabel(frameMiddle,text="Current Missing Value Code:"),
          column=1,row=1,sticky="w")
  tkgrid(tklabel(frameMiddle,text=omiss,width=40),column=2,row=1,sticky="w")
  tkgrid(tklabel(frameMiddle,text="Current nominal level of confidence (p.lev):"),
          column=1,row=2,sticky="w")
  tkgrid(tklabel(frameMiddle,text=oplev,width=40),column=2,row=2,sticky="e")
  tkgrid(tklabel(frameMiddle,text="Segment to which to adjust the series (Iadj):"),
          column=1,row=3,sticky="w")
  tkgrid(tklabel(frameMiddle,text=oadj,width=40),column=2,row=3,sticky="e")
  tkgrid(tklabel(frameMiddle,text="Current Mq (# of points for evaluating PDF):"),
          column=1,row=4,sticky="w")
  tkgrid(tklabel(frameMiddle,text=omq0,width=40),column=2,row=4,sticky="e")
  tkgrid(tklabel(frameMiddle,text="Current Ny4a (max # of years of data for estimating PDF):"),
          column=1,row=5,sticky="w")
  tkgrid(tklabel(frameMiddle,text=oNy4a,width=40),column=2,row=5,sticky="e")
  tkgrid(tklabel(frameMiddle,text="Current input Base series filename:"),
          column=1,row=6,sticky="w")
  tkgrid(tklabel(frameMiddle,text=oact,width=40),column=2,row=6,sticky="w")
  tkgrid(tklabel(frameMiddle,text="Current input Reference series filename:"),
          column=1,row=7,sticky="w")
  tkgrid(tklabel(frameMiddle,text=oref,width=40),column=2,row=7,sticky="w")
  tkgrid(tklabel(frameMiddle,text="Current data directory:    "),
          column=1,row=8,sticky="w")
  tkgrid(tklabel(frameMiddle,text=ocurdir,width=40),column=2,row=8,sticky="e")
  tkgrid(tklabel(frameMiddle,text="Current output directory:    "),
          column=1,row=9,sticky="w")
  tkgrid(tklabel(frameMiddle,text=ooutdir,width=40),column=2,row=9)

  # frameMiddle <- tkframe(frameOverall,relief="groove",borderwidth=2)
  # assign("frameMiddle",frameMiddle,env=.GlobalEnv)
  frameLower <- tkframe(frameOverall,relief="groove",borderwidth=2)
  assign("xscr",tkscrollbar(frameLower,repeatinterval=5,orient="horizontal",
                            command=function(...)tkxview(txt,...)),envir=.GlobalEnv)
  assign("yscr",tkscrollbar(frameLower,repeatinterval=5,
                            command=function(...)tkyview(txt,...)),envir=.GlobalEnv)
  assign("txt",tktext(frameLower,bg="white",font="courier",
                      xscrollcommand=function(...)tkset(xscr,...),
                      yscrollcommand=function(...)tkset(yscr,...)
  ),envir=.GlobalEnv)
  tkgrid(frameUpper)
  tkgrid(frameMiddle)
  tkgrid(txt,yscr)
  tkgrid(xscr)
  tkgrid.configure(yscr,sticky="ns")
  tkgrid.configure(xscr,sticky="ew")
  tkfocus(txt)
  tkgrid(frameLower)
  tkgrid(frameOverall)
}
