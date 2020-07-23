StepSize.wRef<-function(Bseries,Rseries,InCs,output,MissingValueCode,p.lev=0.95,Iadj=10000,Mq=10,GUI=F,Ny4a=0){
  if(Ny4a>0&Ny4a<=5) Ny4a<-5
  if(!p.lev%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
    ErrorMSG<<-paste("FindU: input p.lev",p.lev,"error\n",
                     get("ErrorMSG",env=.GlobalEnv),"\n")
    cat(ErrorMSG)
    return(-1)
  }
  plev<-p.lev
  pkth<-match(p.lev,c(0.75,0.8,0.9,0.95,0.99,0.9999))

  Read.wRef(Bseries,Rseries,MissingValueCode)
  N<-length(Y0); Nadj<-Nt*Ny4a
  # readin PTmax table
  readPTtable(N,pkth)
  # readin Ips
  itmp<-readLines(InCs)
  if(length(itmp)<2) {
    #   stop("There is no input Ips")
    Ns<-0
    Ips<-N
    Ids<-0
  }
  else{
    Ns<-length(itmp)-1
    Ips<-c(rep(0,Ns),N)
    Ids<-rep(0,Ns)
    for(i in 1:Ns){ # using YYYYMMDD as index, searching for the largest
      # date less or equal to given YYYYMMDD
      ymdtmp<-as.numeric(substr(itmp[i+1],7,16))
      it<-match(ymdtmp,IY0)
      if(!is.na(it)) Ips[i]<-it
      else Ips[i]<-max(c(1:N)[IY0<=ymdtmp])
      Ids[i]<-as.numeric(substr(itmp[i+1],1,1))
    }
    if(sum(is.na(Ips))>0|!identical(Ips,sort(Ips))){
      ErrorMSG<<-paste("StepSize.wRef: Ips read in from ",InCs,"error!\n",
                       get("ErrorMSG",env=.GlobalEnv),"\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
    }
  }

  if(Ns>0) {
    Nsegs<-Ips-c(0,Ips[1:Ns])
    Iseg.longest<-sort(Nsegs,index=T,decreasing=T)$ix[1]
  }
  else Iseg.longest<-0

  # if(Iadj==1|Iseg.longest==0) Iseg.adj<-Ns+1 else Iseg.adj<-Iseg.longest
  if(Iadj>(Ns+1)|Iseg.longest==0) Iseg.adj<-Ns+1
  else if(Iadj==0)Iseg.adj<-Iseg.longest
  else Iseg.adj<-Iadj

  otmp<-Rphi(Y0,Ips,Ns)
  W<-otmp$W
  WL<-otmp$WL
  WU<-otmp$WU
  cor<-otmp$cor
  corl<-otmp$corl
  corh<-otmp$corh
  df<-(N-2-Ns)
  p.cor<-pt(abs(cor)*sqrt(df/(1-cor^2)),df)
  ofileIout<-paste(output,"_fCs.txt",sep="")
  ofileMout<-paste(output,"_mCs.txt",sep="")
  ofilePdf<-paste(output,"_F.pdf",sep="")
  ofileSout<-paste(output,"_Fstat.txt",sep="")
  ofileAout<-paste(output,"_F.dat",sep="")
  # ofileRout<-paste(output,"_Base_Ref.fitUDfinal",sep="")
  file.create(ofileSout)
  file.create(ofileAout)
  file.create(ofileIout)
  file.create(ofilePdf)
  # file.create(ofileRout)

  cat(paste("Input Base Series:",Bseries,"\n"),file=ofileSout)
  cat(paste("Input Ref Series:",Rseries,"\n"),file=ofileSout,append=T)
  cat(paste("The adj-diff. autocor is:",round(cor,4),"(",round(corl,4),
            ",",round(corh,4),"p=",round(p.cor,4),")\n"), file=ofileSout,append=T)

  cat(paste(Ns,"changepoints in Series", Bseries,"\n"),
      file=ofileIout)

  if(Ns>0)
    for(i in 1:Ns){
      Ic<-Ips[i]
      Id<-Ids[i]
      I0<-if(i==1) 0 else Ips[i-1]
      I3<-Ips[i+1]
      Nseg<-I3-I0
      PTx95<-getPTx95(cor,Nseg)
      PTx95L<-getPTx95(corl,Nseg)
      PTx95U<-getPTx95(corh,Nseg)

      Pk0<-Pk.PMT(Nseg)
      otmp<-PTKIc(W[(I0+1):I3],Pk0,Ic-I0)
      prob<-otmp$prob
      otmp<-PTKIc(WL[(I0+1):I3],Pk0,Ic-I0)
      probL0<-otmp$prob
      otmp<-PTKIc(WU[(I0+1):I3],Pk0,Ic-I0)
      probU0<-otmp$prob
      probL<-min(probL0,probU0)
      probU<-max(probL0,probU0)

      otmp<-PTKIc(Y0[(I0+1):I3],Pk0,Ic-I0)
      PTx0<-otmp$PTk
      if(Id==0) { # type-0 changepoints
        if(probU<plev) Idc<-"No  "
        else if(probL<plev&probU>=plev) Idc<-"?   "
        else if(probL>=plev) Idc<-"YifD"
        if(PTx0>=PTx95U) Idc<-"Yes "
      }
      else if(Id==1) { # type-1 changepoints
        if(PTx0<PTx95L) Idc<-"No  "
        else if(PTx0>=PTx95L&PTx0<PTx95U) Idc<-"?   "
        else if(PTx0>=PTx95U) Idc<-"Yes "
      }
      cat(paste(sprintf("%1.0f",Id)," ",
                sprintf("%-4.4s",Idc),
                sprintf("%10.0f",IY0[Ic])," (",
                sprintf("%10.4f",probL),"-",
                sprintf("%10.4f",probU),")",
                sprintf("%6.3f",plev),
                sprintf("%10.4f",PTx0)," (",
                sprintf("%10.4f",PTx95L),"-",
                sprintf("%10.4f",PTx95U),")\n",sep=""),
          file=ofileIout,
          append=TRUE)
      cat(paste("PMT : c=", sprintf("%4.0f",Ic),
                "; (Time ", sprintf("%10.0f",IY0[Ic]),
                "); Type= ",sprintf("%4.0f",Id),
                "; p=", sprintf("%10.4f",prob),
                "(", sprintf("%10.4f",probL),
                "-", sprintf("%10.4f",probU),
                "); PTmax=", sprintf("%10.4f",PTx0),
                "; CV95=", sprintf("%10.4f",PTx95),
                "(", sprintf("%10.4f",PTx95L),
                "-", sprintf("%10.4f",PTx95U),
                "); Nseg=", sprintf("%4.0f",Nseg),
                "\n",sep=""), file=ofileSout, append=T)
    }

  # estimate delta from Y0 (Base-Ref)
  otmp<-Rphi(Y0,Ips,Ns)
  cor<-otmp$cor
  muDif<-rep(0,Ns+1)
  Ro<-Y0
  Wo<-Y0
  omuDif<-Y0
  oY0<-Y0
  for(i in 1:(Ns+1)){
    I0<- if(i==1) 1 else Ips[i-1]+1
    I2<- if(i>Ns) length(Y0) else Ips[i]
    muDif[i]<-mean(Y0[I0:I2])
    omuDif[I0:I2]<-muDif[i]
    Ro[I0:I2]<-Y0[I0:I2]-muDif[i]
  }
  Wo[1]<-Ro[1]
  Wo[2:N]<-Ro[2:N]-cor*Ro[1:(N-1)]*IY0flg[1:(N-1)]
  # write.table(cbind(IY0,round(oY0,4),round(omuDif,4),round(Wo,4)),
  #             file=ofileRout,col.names=F,row.names=F)

  # transfer Ips(Base-Ref) to Ips(Base)
  Ips0<-Ips
  IY1<-bdata[,1]*10000+bdata[,2]*100+bdata[,3]
  IYM<-bdata[,2]*100+bdata[,3]
  assign('Imd',IYM,env=.GlobalEnv)
  inds<-sort(unique(IYM))
  rtmp<-cbind(1:length(IY1),IY1)
  Ips<-merge(IY0[Ips0],rtmp,by.x=1,by.y="IY1")[,2]
  Ips[Ns+1]<-length(IY1)

  Ti<-TiB
  assign("Ti",Ti,envir=.GlobalEnv)
  N<-length(TiB)
  IY0flg<-IYBflg
  assign("IY0flg",IY0flg,envir=.GlobalEnv)

  dtmp<-rmCycle(bdata)
  adjBase<-dtmp$Base
  EB<-rep(0,length(IY1))
  for(i in 1:length(IY1))
    EB[i]<-dtmp$EB[inds==IYM[i]]

  if(Ns>0)
    for(i in 1:(Ns+1)){
      I0<- if(i==1) 1 else Ips[i-1]+1
      I2<- if(i>Ns) N else Ips[i]
      DeltaD<-muDif[i]-muDif[Iseg.adj]
      adjBase[I0:I2]<-adjBase[I0:I2]+EB[I0:I2]-DeltaD # adjBase is base series adj to Iseg.adj
    }
  EPBa<-mean(adjBase)

  Ydata<-cbind(bdata[,1:3],adjBase)
  dtmp<-rmCycle(Ydata)  # adjusted and de-seasonalized Base series
  EB0<-dtmp$EB
  EEBd<-mean(EB0)
  EB<-rep(0,length(IY1))
  for(i in 1:length(IY1))
    EB[i]<-dtmp$EB[inds==IYM[i]]
  Aadj<-dtmp$Base  # de-seasonalize Badj
  Ipd<-c(N)
  # dtmp<-LSmultipleRed(Aadj,Ti,Ipd)
  assign('EB',EB0,env=.GlobalEnv)
  dtmp<-LSmultiRedCycle(Aadj,Ti,Ipd,Iseg.adj)
  muD<-dtmp$mu[1]
  betaD<-dtmp$trend
  betaDL<-dtmp$betaL
  betaDU<-dtmp$betaU
  corD<-dtmp$cor
  corDL<-dtmp$corl
  corDU<-dtmp$corh
  p.trD<-dtmp$p.tr

  dtmp<-rmCycle(bdata) # de-seasonalized Base series
  assign('EB',dtmp$EB,env=.GlobalEnv)
  tbase<-dtmp$Base
  Ipd<-length(tbase)
  # dtmp<-LSmultipleRed(tbase,Ti,Ipd)
  dtmp<-LSmultiRedCycle(tbase,Ti,Ipd,Iseg.adj)
  beta0<-dtmp$trend
  meanhat0<-dtmp$meanhat
  Ehat0<-mean(meanhat0)
  cor<-dtmp$cor
  corL<-dtmp$corl
  corU<-dtmp$corh

  cat(paste("Ignore changepoints -> trend0 =",round(beta0,6),
            "(",round(dtmp$betaL,6),",",round(dtmp$betaU,6),
            ") (p=",round(dtmp$p.tr,4),
            "); cor=", round(cor,4),"(",round(corL,4),",",
            round(corU,4),")\n\n"),
      file=ofileSout,append=TRUE)
  cat("Step-sizes estimated from difference series:\n",
      file=ofileSout,append=TRUE)
  if(Ns>0) cat(round(muDif[2:(Ns+1)]-muDif[1:Ns],4),
               file=ofileSout,append=TRUE,fill=80)
  cat(paste("\n after such adjustments, the base series trend=",
            round(betaD,6),"(",round(betaDL,6),",",round(betaDU,6),
            ") (p=",round(p.trD,4),"); cor=",round(corD,3),"(",round(corDL,3),
            ",",round(corDU,3),")\n\n"),file=ofileSout,append=TRUE)

  otmp<-rmCycle(bdata)
  EB<-rep(0,length(IY1))
  for(i in 1:length(IY1))
    EB[i]<-otmp$EB[inds==IYM[i]]
  Base<-otmp$Base
  EPB<-mean(bdata[,"data"],na.rm=T)
  tt<-TRUE
  Niter<-0
  while(tt){
    Niter<-Niter+1
    EB00<-EB
    otmp<-LSmultipleRed(Base,Ti,Ips)
    meanhat<-otmp$meanhat
    df<-(N-2-Nt-Ns)
    p.cor<-pt(abs(otmp$cor)*sqrt(df/(1-otmp$cor^2)),df)
    otrend<-otmp$trend
    corout<-c(otmp$cor,otmp$corl,otmp$corh,p.cor)
    sig<-otmp$sig
    p.tro<-otmp$p.tr
    if(Ns==0) {
      EB1<-EB00
      tt<-FALSE
    }
    else{
      for(i in 1:(Ns+1)){
        I0<- if(i==1) 1 else Ips[i-1]+1
        I2<- if(i>Ns) length(Base) else Ips[i]
        Delta<- sig[i+1]-sig[Iseg.adj+1]
        Base[I0:I2]<-Base[I0:I2]+EB00[I0:I2]-Delta
      }
      # re-estimate seasonal cycle using adjusted series:
      EB1<-rep(0,length(inds))
      for(i in 1:length(inds))
        EB1[i]<-mean(Base[IYM==inds[i]],na.rm=T)
      for(i in 1:length(IY1))
        EB[i]<-EB1[inds==IYM[i]]
      for(i in 1:(Ns+1)){
        I0<- if(i==1) 1 else Ips[i-1]+1
        I2<- if(i>Ns) length(Base) else Ips[i]
        Delta<- sig[i+1]-sig[Iseg.adj+1]
        Base[I0:I2]<-Base[I0:I2]-EB[I0:I2]+Delta
      }
      VEB<-sqrt(var(EB1))
      if(is.na(VEB)) tt<-FALSE
      else { if(max(abs(EB0-EB1))<VEB/1000|Niter>20) tt<-FALSE }
      EB0<-EB1
    }
  }
  Ehat<-mean(meanhat)
  meanhat0<-meanhat0-Ehat0+Ehat
  Ro<-Base-meanhat
  Ro[2:N]<-Ro[2:N]-corout[1]*Ro[1:(N-1)]

  IY1<-bdata[,1]*10000+bdata[,2]*100+bdata[,3]

  if(Ns>0){
    #   Rb<-Base-otmp$trend*Ti+EB
    #   QMout<-QMadjGaussian(Rb,Ips,Mq,Iseg.adj,Nadj)
    QMout<-QMadjGaussian.wRef(bdata[,4],bdata[,4]-bdata[,5],Ips,Mq,Iseg.adj,Nadj,Nt,Ny4a)
    B<-QMout$PA
    cat(paste("Nseg_shortest =",QMout$Nseg.mn,"; Mq = ",QMout$Mq,"; Ny4a = ",Ny4a,"\n"),
        file=ofileSout,append=T)
    cat(paste("\n Adjust to segment", Iseg.adj,": from",
              if(Iseg.adj==1) 1 else Ips[Iseg.adj-1]+1,
              "to",Ips[Iseg.adj],"\n"),file=ofileSout,append=T)
    #   cat("#Fcat, DP (CDF and Differnces in category mean)\n",file=ofileSout,
    #       append=T)
    if(QMout$Mq>1){
      oline<-paste('#Fcat: frequency category boundaries\n',
                   '#DP: Difference in the category means\n#',sep='')
      for(i in 1:Ns) oline<-paste(oline,paste('Fcat.Seg',i,sep=''),paste('DP.Seg',i,sep=''))
      oline<-paste(oline,'\n')
      cat(oline,file=ofileSout,append=T)

      write.table(round(QMout$osmean,4),file=ofileSout,append=T,
                  row.names=F,col.names=F)
      for(i in 1:(Ns+1)){
        I1<-if(i==1) 1 else Ips[i-1]+1
        I2<-Ips[i]
        if(i!=Iseg.adj)
          cat(paste("Seg. ",i,": mean of QM-adjustments =",round(QMout$AdjM[i],4),
                    "\n",sep=""),file=ofileSout,append=T)
      }
    }
  }
  else B<-Base
  # else B<-Base-otmp$trend*Ti+EB
  # B<-B+otmp$trend*Ti

  pdf(file=ofilePdf,onefile=T,paper='letter')
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  par(mfrow=c(3,1))
  par(mar=c(3,4,3,2)+.1,mgp=c(1,.5,0),cex.main=.8,cex.lab=.8,cex.axis=.8,cex=.8)

  uyrs<-unique(floor(ori.bdata[,1]/10))*10
  labels<-NULL
  ats<-NULL
  for(i in 1:length(uyrs)){
    if(!is.na(match(uyrs[i],ori.bdata[,1]))){
      labels<-c(labels,uyrs[i])
      ats<-c(ats,match(uyrs[i],ori.bdata[,1]))
    }
  }

  adj<-Base+EB
  adjB<-Base+EB
  meanhatD<-rep(0,N)
  if(Ns>0) for(i in 1:(Ns+1)){
    I1<-if(i==1) 1 else Ips[i-1]+1
    I2<-Ips[i]
    Delta<-sig[Iseg.adj+1]-sig[i+1]
    DeltaD<-muDif[Iseg.adj]-muDif[i]
    adj[I1:I2]<-adj[I1:I2]+Delta
    adjB[I1:I2]<-adjB[I1:I2]+DeltaD
    meanhatD[I1:I2]<-EEBd+muD+betaD*Ti[I1:I2]-DeltaD
  }
  else meanhatD<-EPB+muD+betaD*Ti

  IYori<-ori.bdata[,1]*10000+ori.bdata[,2]*100+ori.bdata[,3]
  rtmp<-cbind(IY0,oY0,omuDif)
  stmp<-merge(rtmp,t(t(IYori)),all.y=TRUE,by.x="IY0",by.y=1)
  pdata<-stmp[,2]

  plot(1:dim(ori.bdata)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(oY0),max(oY0)),
       xaxt="n",col="black",lwd=0.5,
       main="a. Base-minus-reference series")
  axis(side=1,at=ats,labels=labels)

  pdata<-stmp[,3]
  lines(1:dim(ori.bdata)[1],pdata,col="red")

  pdata[owflg]<-Base
  plot(1:dim(ori.bdata)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(Base,otmp$meanhat),max(Base,otmp$meanhat)),
       xaxt="n",col="black",lwd=0.5,
       main="b. De-seasonalized base series")
  axis(side=1,at=ats,labels=labels)
  pdata[owflg]<-otmp$meanhat
  lines(1:dim(ori.bdata)[1],pdata,col="red")

  pdata[owflg]<-Base+EB
  plot(1:dim(ori.bdata)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(Base+EB),max(Base+EB)),
       xaxt="n",col="black",lwd=0.5,
       main="c. Base series")
  axis(side=1,at=ats,labels=labels)
  pdata[owflg]<-otmp$meanhat+mean(EB1)
  lines(1:dim(ori.bdata)[1],pdata,col="red")
  pdata[owflg]<-meanhatD
  lines(1:dim(ori.bdata)[1],pdata,col="blue")

  xr<-dim(ori.bdata)[1]*.1; yrr<-mean((Base+EB)[1:xr])
  if(abs(yrr-max(Base+EB))<abs(yrr-min(Base+EB))) {
    yrr1<-min(Base+EB)+.2*(max(Base+EB)-min(Base+EB))
    yrr2<-min(Base+EB)+.05*(max(Base+EB)-min(Base+EB))
  }
  else {
    yrr1<-max(Base+EB)-.05*(max(Base+EB)-min(Base+EB))
    yrr2<-max(Base+EB)-.2*(max(Base+EB)-min(Base+EB))
  }
  text(2,yrr1,"Adjustments estimated from Base",adj=0,cex=.8); lines(dim(ori.bdata)[1]*c(.31,.33),rep(yrr1,2),col="red")
  text(2,yrr2,"Adjustments estimated from Diff",adj=0,cex=.8); lines(dim(ori.bdata)[1]*c(.31,.33),rep(yrr2,2),col="blue")

  pdata[owflg]<-adj
  plot(1:dim(ori.bdata)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(c(adj,B)),max(c(adj,B))),
       xaxt="n",col="black",lwd=0.5,
       main="d. Mean-adjusted base series")
  axis(side=1,at=ats,labels=labels)
  pdata[owflg]<-adjB
  lines(1:dim(ori.bdata)[1],pdata,col="red")

  xr<-dim(ori.bdata)[1]*.1; yrr<-mean(adj[1:xr])
  if(abs(yrr-max(c(adj,B)))<abs(yrr-min(c(adj,B)))) {
    yrr1<-min(c(adj,B))+.2*(max(c(adj,B))-min(c(adj,B)))
    yrr2<-min(c(adj,B))+.05*(max(c(adj,B))-min(c(adj,B)))
  }
  else {
    yrr1<-max(c(adj,B))-.05*(max(c(adj,B))-min(c(adj,B)))
    yrr2<-max(c(adj,B))-.2*(max(c(adj,B))-min(c(adj,B)))
  }
  text(2,yrr1,"Adjustments estimated from Base",adj=0,cex=.8); lines(dim(ori.bdata)[1]*c(.31,.33),rep(yrr1,2),col="black")
  text(2,yrr2,"Adjustments estimated from Diff",adj=0,cex=.8); lines(dim(ori.bdata)[1]*c(.31,.33),rep(yrr2,2),col="red")

  if(Ns>0) if(QMout$Mq>1){
    pdata[owflg]<-B
    plot(1:dim(ori.bdata)[1],pdata,type="l",xlab="",ylab="",
         ylim=c(min(c(adj,B)),max(c(adj,B))),
         xaxt="n",col="black",lwd=0.5,
         main="e. QM-adjusted base series")
    axis(side=1,at=ats,labels=labels)
  }

  # test plot
  if(Ns>0) if(QMout$Mq>1){
    par(mar=c(4,5,3,2)+.1,cex=.8,mfrow=c(2,2),mgp=c(1.2,.5,0))
    col=0
    np<-0
    osp<-QMout$osp
    osmean<-QMout$osmean
    for(i in 1:(Ns+1)){
      Fd<-.5/QMout$Mq
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i]
      ymax<-max(osp[,2:3],na.rm=T); ymin<-min(osp[,2:3],na.rm=T)
      if(i!=Iseg.adj){
        np<-np+1
        if(col==0) {
          col<-2
          plot(osp[I1:I2,2],osp[I1:I2,3],xlim=c(0,1),ylim=c(ymin,ymax),
               type="l",lwd=1,col=col,xlab="Cumulative Frequency",
               ylab="QM Adjustment")
          title(cex.main=.9,main=paste("distribution of QM adjustments with Mq=",QMout$Mq),line=.5)
          icol<-2*np
          for(j in 1:QMout$Mq){
            lines(c(osmean[(j+1),icol]-Fd,osmean[(j+1),icol]+Fd),
                  c(rep(osmean[(j+1),(icol+1)],2)),col=col,lty=2,lwd=.5)
            if(j>=1&j<QMout$Mq) lines(rep(osmean[(j+1),icol]+Fd,2),
                                      c(osmean[(j+1),(icol+1)],osmean[(j+2),(icol+1)]),col=col,lty=2,lwd=.5)
          }
        }
        else{
          col<-col+1
          lines(osp[I1:I2,2],osp[I1:I2,3],lwd=1,col=col)
          icol<-2*np
          for(j in 1:QMout$Mq){
            lines(c(osmean[(j+1),icol]-Fd,osmean[(j+1),icol]+Fd),
                  c(rep(osmean[(j+1),(icol+1)],2)),col=col,lty=2,lwd=.5)
            if(j>=1&j<QMout$Mq) lines(rep(osmean[(j+1),icol]+Fd,2),
                                      c(osmean[(j+1),(icol+1)],osmean[(j+2),(icol+1)]),col=col,lty=2,lwd=.5)
          }
        }
        text(.15,ymax-np*(ymax-ymin)/(Ns*3),paste("Seg.",i))
        lines(c(.25,.30),rep(ymax-np*(ymax-ymin)/(Ns*3),2),lwd=2,col=col)
      }
      else np<-np+1
    }
  }

  par(op)
  dev.off()
  cat("Common trend TPR fit to the de-seasonalized Base series:\n",
      file=ofileSout,append=TRUE)
  cat(paste("#steps= ",Ns,"; trend=",round(otrend,6),"(p=",
            round(p.tro,4),"); cor=",
            round(corout[1],4),"(",round(corout[2],4),",",round(corout[3],4),
            ")  p=",round(corout[4],4), "\n"),
      file=ofileSout,append=TRUE)
  oout<-NULL
  if(Ns>0){
    for(i in 1:Ns){
      stepsize<-if(i==1) sig[i+2] else sig[i+2]-sig[i+1]
      oout<-c(oout,stepsize)
    }
    cat(round(oout,4),file=ofileSout,append=TRUE,fill=80)
  }

  odata<-matrix(NA,dim(ori.bdata)[1],12)
  odata[owflg,1]<-Ti
  odata[,2]<-ori.bdata[,1]*10000+ori.bdata[,2]*100+ori.bdata[,3]
  # odata[owflg,3]<-round(Base+EB,4)
  odata[,3]<-ori.bdata[,4]
  odata[owflg,4]<-round(meanhatD,4)
  odata[owflg,5]<-round(adjB,4)
  odata[owflg,6]<-round(otmp$meanhat+mean(EB1),4)
  odata[owflg,7]<-round(adj,4)
  odata[owflg,8]<-round(Base,4)
  odata[owflg,9]<-round(otmp$meanhat,4)
  odata[owflg,10]<-round(otmp$meanhat+EB,4)
  # odata[owflg,11]<-round(Ro,4)
  if(Ns>0) if(QMout$Mq>1) odata[owflg,11]<-round(B,4)
  odata[owflg,12]<-round(meanhat0,4)

  Imd1<-ori.bdata[,2]*100+ori.bdata[,3]
  if(sum(is.na(ori.bdata[,4])==F&Imd1==229)>0){
    if(Ns>0){
      tdata<-ori.bdata[is.na(ori.bdata[,4])==F,]
      IY1<-tdata[,1]*10000+tdata[,2]*100+tdata[,3]
      Ips.ymd<-IY0[Ips0]
      Ips.1<-rep(NA,Ns+1)
      for(i in 1:Ns) Ips.1[i]<-c(1:length(IY1))[IY1==Ips.ymd[i]]
      Ips.1[Ns+1]<-length(IY1)
      #     Ips.1<-c(1:length(IY1))[Ips.ymd==IY1]
      Imd2<-tdata[,2]*100+tdata[,3]
      Ids.leap<-c(1:length(Imd2))[Imd2==229]
      Nl<-length(Ids.leap)
      Rb<-Base-otmp$trend*Ti+EB
      Rb1<-tdata[,4]; Rb1[-Ids.leap]<-Rb
      Ti1<-rep(NA,length(IY1)); Ti1[-Ids.leap]<-Ti
      for(i in 1:length(Ids.leap)) {
        Rb1[Ids.leap[i]]<-tdata[Ids.leap[i],4]+Rb1[Ids.leap[i]-1]-tdata[Ids.leap[i]-1,4]
        Ti1[Ids.leap[i]]<-Ti1[Ids.leap[i]-1]
      }
      if(QMout$Mq>1){
        B1<-QMadjGaussian(Rb1,Ips.1,Mq,Iseg.adj,Nadj)$PA
        B1<-B1+otmp$trend*Ti1
        B1.leap<-B1[Ids.leap]
        odata[is.na(odata[,3])==F&Imd1==229,11]<-round(B1.leap,4)
      }
    }
    else
      odata[Imd1==229,11]<-odata[Imd1==229,3]
    Ids.leapo<-c(1:dim(ori.bdata)[1])[is.na(ori.bdata[,4])==F&Imd1==229]
    for(jth in 1:length(Ids.leapo)){
      kth<-Ids.leapo[jth]
      if(Ns>0){
        k1th<-if(odata[kth-1,2]%in%IY0[Ips0]) (kth+1) else (kth-1)
      }
      else k1th<-kth-1
      for(pth in c(4,6,9,10,12)) odata[kth,pth]<-odata[k1th,pth]
      for(pth in c(5,7,8)){delta1<-odata[k1th,3]-odata[k1th,pth]; odata[kth,pth]<-odata[kth,3]-delta1}
    }
  }

  write.table(file=ofileAout,odata,col.names=F,row.names=F,na=MissingValueCode)
  if(GUI) return(0)
  else {
    file.copy(from=ofileIout,to=ofileMout,overwrite=TRUE)
    cat("StepSize.wRef finished successfully...\n")
  }
}
