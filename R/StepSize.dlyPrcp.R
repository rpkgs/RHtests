StepSize.dlyPrcp<-function(InSeries,InCs,output,
  MissingValueCode="-999.99",GUI=FALSE, pthr=0,Mq=10,p.lev=0.95,Iadj=10000,Ny4a=0)
{
  ErrorMSG<-NA
  Ncat.min<-20
  if(Ny4a>0&Ny4a<=5) Ny4a<-5
  assign("ErrorMSG",ErrorMSG,envir=.GlobalEnv)
  if(!p.lev%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
      ErrorMSG<<-paste("StepSize.dlyprcp: input p.lev",p.lev,"error\n",
                 get("ErrorMSG",env=.GlobalEnv),"\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
  }
  plev<-p.lev
  pkth<-match(p.lev,c(0.75,0.8,0.9,0.95,0.99,0.9999))
  itmp<-ReadDLY(InSeries,MissingValueCode,pthr)
  N<-length(Y0)
  Nall<-dim(ori.itable)[1]; Pfreq<-N/Nall; Nadj<-Ny4a*366*Pfreq
  readPFtable(N,pkth)
  if(itmp==(-1)){
    ErrorMSG<<-paste("StepSize.dlyprcp: Error in read data from",InSeries,"\n",
               get("ErrorMSG",env=.GlobalEnv),"\n")
    if(!GUI) cat(ErrorMSG)
    return(-1)
  }
  itmp<-readLines(InCs)
  if(length(itmp)<2){
    Ns<-0
    Ips<-N
  }  else{
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
      ErrorMSG<<-paste("StepSize.dlyprcp: Ips read in from ",InCs,"error:")
      for(i in 1:Ns)
        ErrorMSG<<-paste(get("ErrorMSG",env=.GlobalEnv),Ips[i])
      ErrorMSG<<-paste(get("ErrorMSG",env=.GlobalEnv),"\n\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
    }
  }
# smallestP<-min(Y0)
  smallestP<-min(ori.itable[ori.itable[,"data"]>0,"data"],na.rm=T)
  smallestP<-max(c(smallestP,pthr))
  sumLogPi<-sum(log(Y0))

  ofileSout<-paste(output,"_Fstat.txt",sep="")
  file.create(ofileSout)

  cat(paste("Input data filename:","; N=",N,"; Nday=",Nall,
      "; pthr=",pthr,"\n"), file=ofileSout)
  cat(paste("The nominal level of confidence (1-alpha)=",plev,"\n"),
      file=ofileSout,append=T)
  cat(paste("The smallest non-zero dailyP is:",smallestP,"\n"),
      file=ofileSout,append=T)

  Nsegs<-Ips-c(0,Ips[1:Ns])
  Iseg.longest<-sort(Nsegs,index=T,decreasing=T)$ix[1]
  if(Iadj>(Ns+1)|Iseg.longest==0) Iseg.adj<-Ns+1
  else if(Iadj==0)Iseg.adj<-Iseg.longest
  else Iseg.adj<-Iadj

  P<-Y0
  aa<-(-1)
  bb<-1
  lambdas<-seq(-1,1,.01)
  LHs<-rep(NA,length(lambdas))
  Niter0<-0
  LH.max<-(-9999999.)
  lr.max<-(-1)
  for(dr in c(1,0.5,0.2,0.1)){
    lrs<-seq(aa,bb,by=dr)
    for(lr in 1:length(lrs)){
      lambda<-aa+(lr-1)*dr
      ind<-round((lambda+1)*100+1)
      if(abs(lambdas[ind]-lambda)>1e-8) {
        print(c(lambdas[ind],lambda,lambdas[ind]-lambda))
        stop(paste("ind=",ind,"lambda=",lambda,"error"))
      }
      if(is.na(LHs[ind])==F){  # already done before
        LH<-LHs[ind]
        Y1<-BCtrans(P,lambda)
      }
      else{
        Y1<-BCtrans(P,lambda)
	otmp<-LSmultipleRed(Y1,Ti,Ips)
	SSEf<-sum(otmp$resi**2)
	LH<--log(SSEf*(N-3)/(N-2-Ns))*N/2+(lambda-1)*sumLogPi
        LHs[ind]<-LH
      }
      if(LH>LH.max){
        LH.max<-LH
	lambda.max<-lambda
      }
    }
    if(dr>0.1){
      aa<-max(c(lambda.max-dr,-1))
      bb<-min(c(lambda.max+dr,1))
      lrs<-seq(aa,bb,by=dr)
    }
  }
  lambda<-lambda.max
  Y1<-BCtrans(Y0,lambda)

  Ipd<-c(N)
  dtmp<-LSmultipleRed(Y1,Ti,Ipd)
  beta0<-dtmp$trend
  cat(paste("Ignore changepoints -> trend0 =",round(beta0,6),
    "(p=", round(dtmp$p.tr,4), "); cor=", round(dtmp$cor,4),
    "(",round(dtmp$corl,4),",",round(dtmp$corh,4),")\n\n"),
    file=ofileSout,append=TRUE)
  
  otmp<-LSmultipleRed(Y1,Ti,Ips)
  cor<-otmp$cor
  corL<-otmp$corl
  corU<-otmp$corh
  W<-otmp$W
  WL<-otmp$WL
  WU<-otmp$WU

  Rresi<-otmp$resi
  pcor<-pt(abs(cor)*sqrt((N-2)/(1-cor^2)),N-2)
  meanhat<-otmp$meanhat
  meanhatA<-rep(NA,N)
  Pmeanhat<-IVBCtrans(meanhat,lambda)

# cat("Common trend TPR fit to the transformed Base series:\n",
#     file=ofileSout,append=T)
# cat(paste("#steps= ",Ns,"; trend=",round(otmp$trend,6),"unit/yr (p=",
#           round(otmp$p.tr,4),"); lambda=",round(lambda,2),"; cor=",
#           round(cor,4),"(",round(corL,4),",",round(corU,4),") p=",
#           round(pcor,4),"\n"),file=ofileSout,append=T)
  if(Ns>0) {
    cat(paste("Step-sizes of the transfromed/original series:\n"),file=ofileSout,append=T)
    C<-rep(NA,Ns); sumC<-0; E<-rep(0,Ns+1)
    for(i in 1:Ns){
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i+1]
      Ic<-Ips[i]
      meanhat0<-otmp$mu[i]+otmp$trend*Ti[I1:I2]
      Y<-otmp$mu[i+1]+otmp$trend*Ti[I1:I2]
      R0<-IVBCtrans(meanhat0,lambda)
      R<-IVBCtrans(Y,lambda)
      C[i]<-R[Ic-I1+1]-R0[Ic-I1+1]
      sumC<-sumC+C[i]
      E[i+1]<-sumC
      stepsize<-otmp$mu[i+1]-otmp$mu[i]
      cat(paste(Ips[i],IY0[Ips[i]],
                "transfered stepsize=",round(stepsize,4),
		"original stepsize=",round(C[i],4),
		"\n"), file=ofileSout,append=T)
    }
# calculate IBC adjustment
    Y<-rep(NA,N)
    for(i in 1:(Ns+1)){
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i]
      Delta<-otmp$mu[Iseg.adj]-otmp$mu[i]
      Y[I1:I2]<-Y1[I1:I2]+Delta
      meanhatA[I1:I2]<-meanhat[I1:I2]+Delta
    }
#   PmeanhatA<-IVBCtrans(meanhatA,lambda)
    PmeanhatA<-rep(NA,N)
    dPmu<-PmeanhatA[1]-Pmeanhat[1]
    PdA<-rep(NA,N)
    for(i in 1:(Ns+1)){
      I1<-if(i==1) 1 else Ips[i-1]+1
      I2<-Ips[i]
      if(i==Iseg.adj){
        PdA[I1:I2]<-P[I1:I2]
        DeltaP<-0
	PmeanhatA[I1:I2]<-Pmeanhat[I1:I2]
      }
      else{
        DeltaP<-E[Iseg.adj]-E[i]
	DeltaP0<-DeltaP
        PdA[I1:I2]<-P[I1:I2]+DeltaP
        PmeanhatA[I1:I2]<-Pmeanhat[I1:I2]+DeltaP
	PdA[I1:I2][PdA[I1:I2]<smallestP]<-smallestP
	Diff<-mean(PdA[I1:I2]-P[I1:I2])
	Delta<-DeltaP-Diff
	tflg<-(Delta<=(-0.01))
	Niter<-1
	while(tflg){
	  Niter<-Niter+1
	  Delta0<-Delta
	  PdA[I1:I2]<-PdA[I1:I2]+Delta
	  PmeanhatA[I1:I2]<-PmeanhatA[I1:I2]+Delta
	  PdA[I1:I2][PdA[I1:I2]<smallestP]<-smallestP
	  Diff<-mean(PdA[I1:I2]-P[I1:I2])
	  Delta<-DeltaP-Diff
	  tflg<-(Delta<=(-0.01)&abs(Delta)<abs(Delta0)&Niter<5)
	}
      }
    }
#   PdA[PdA<smallestP]<-smallestP
    Ptr.mx<-max(PmeanhatA)

    Pdtr<-P+Ptr.mx-PmeanhatA
    
    QMout<-adjDLYp(Pdtr,Ips,Mq,Iseg.adj,Ptr.mx,PmeanhatA,Ncat.min,Nadj)
#   assign("QMout",QMout,env=.GlobalEnv)
    PA<-QMout$PA
#   PA[PA<smallestP]<-max(c(smallestP,pthr))
    PA[PA<smallestP]<-smallestP
    cat(paste("Nseg_shortest =",QMout$Nseg.mn,"; Mq = ",QMout$Mq,"\n"),
        file=ofileSout,append=T)
    cat(paste("\n Adjust to segment", Iseg.adj,": from",
        if(Iseg.adj==1) 1 else Ips[Iseg.adj-1]+1,
        "to",Ips[Iseg.adj],"\n"),file=ofileSout,append=T)

    if(QMout$Mq>1){
      cat("#Fcat, DP (CDF and Differnces in category mean)\n",file=ofileSout,
          append=T)
      write.table(round(QMout$osmean,4),file=ofileSout,append=T,
                  row.names=F,col.names=F)
      for(i in 1:(Ns+1)){
        I1<-if(i==1) 1 else Ips[i-1]+1
        I2<-Ips[i]
        if(i!=Iseg.adj)
        cat(paste("Seg. ",i,": mean of QM-adjustments =",round(mean(QMout$W[I1:I2]),4),
            "\n",sep=""),file=ofileSout,append=T)
      }
    }
    if(QMout$Mq==1) PA<-rep(NA,length(QMout$PA))

    odata<-cbind(1:nrow(ori.itable),ori.itable[,1]*10000+ori.itable[,2]*100+ori.itable[,3],
#                ori.itable[,4],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
                 ori.itable[,4],ori.itable[,4],ori.itable[,4],ori.itable[,4],ori.itable[,4],ori.itable[,4],
                 ori.itable[,4],ori.itable[,4],ori.itable[,4],0,0)
    Imd<-itable[,1]
    odata[Imd,4]<-round(Pmeanhat,4)
    odata[Imd,5]<-round(PA,4)
    odata[Imd,6]<-round(PdA,4)
    odata[Imd,7]<-round(PmeanhatA,4)
    odata[Imd,8]<-round(Y1,4)
    odata[Imd,9]<-round(meanhat,4)
    odata[Imd,10]<-round(Y,4)
    odata[Imd,11]<-round(meanhatA,4)
    odata[Imd,12]<-round(PA-itable[,5],4)
    odata[Imd,13]<-round(PdA-itable[,5],4)
    odataP<-cbind(itable[,1],IY0,itable[,5],Pmeanhat,PA,PdA,PmeanhatA,
                  Y1,meanhat,Y,meanhatA,PA-itable[,5],PdA-itable[,5])
    otrend.ori<-getMtrendFdly(itable[,2:5])
    otrend.QM<-getMtrendFdly(cbind(itable[,2:4],PA))
    otrend.IBC<-getMtrendFdly(cbind(itable[,2:4],PdA))
    cat(paste("Linear trend in the original monthly total P is trend0=",
        round(otrend.ori,5),"mm/month\n"),file=ofileSout,append=T)
    cat(paste("Linear trend in QMadjusted monthly total P is trendQM=",
        round(otrend.QM,5),"mm/month\n"),file=ofileSout,append=T)
    cat(paste("Linear trend in mean adjusted(IBC) monthly total P is trend(IBC)=",
        round(otrend.IBC,5),"mm/month\n\n"),file=ofileSout,append=T)
  }
  else{
    odataP<-cbind(itable[,1],IY0,itable[,5],Pmeanhat,itable[,5],itable[,5],
                  Pmeanhat,Y1,meanhat,Y1,meanhat,0,0)
    odata<-cbind(1:nrow(ori.itable),ori.itable[,1]*10000+ori.itable[,2]*100+ori.itable[,3],
                 ori.itable[,4],ori.itable[,4],0,0,0,0,0,0,0,0)
    otrend.ori<-getMtrendFdly(itable[,2:5])
    cat(paste("Linear trend in the original monthly total P is trend0=",
        round(otrend.ori,5),"mm/month\n"),file=ofileSout,append=T)
  }
# ofileDout<-paste(output,"_adj_F.dat",sep="")
# write.table(odata,file=ofileDout,na=MissingValueCode,col.names=F,row.names=F)

  ofilePout<-paste(output,"_F.dat",sep="")
# write.table(round(odataP,4),file=ofilePout,na=MissingValueCode,col.names=F,row.names=F)
  write.table(odata,file=ofilePout,na=MissingValueCode,col.names=F,row.names=F)

  ofilePdf<-paste(output,"_F.pdf",sep="")
  pdf(file=ofilePdf,onefile=T)
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  par(mfrow=c(2,1))
# par(mar=c(3,4,3,2)+.1)
  par(mar=c(2,4,3,1)+.1)
  par(cex=.8,cex.main=.8,cex.lab=.8,cex.axis=.8)

  if(Ns>0){
    p1data<-cbind(c(1:Nall),ori.itable[,4],NA)
    p1data[Imd,2]<-P
    p1data[Imd,3]<-odataP[,4]

    for(Iseg in 1:Ns){ # plot P~Pmeanhat
      I1<-max(c(1,(itable[Ips[Iseg],1]-180)))
      Ic<-itable[Ips[Iseg],1]
      I2<-min(c(Nall,(itable[Ips[Iseg],1]+180)))
      plot(I1:I2,ori.itable[I1:I2,4],type="l",col="black",ylab="prcp(mm)",
           main=paste("Original dailyP>pthr series\n","Changepoint ",Iseg,
	   " at:",IY0[Ips[Iseg]],sep=""),xaxt="n")
      ats<-c(I1,Ic,I2)
      labels=ori.itable[ats,1]*10000+ori.itable[ats,2]*100+ori.itable[ats,3]
      axis(side=1,at=c(I1,Ic,I2),labels=labels)

      tflg<-is.na(p1data[I1:I2,3])==F
      lines(c(I1:I2)[tflg],p1data[I1:I2,3][tflg],col="red",lwd=1.5)
    }
    for(Iseg in 1:Ns){ # plot Y1~meanhat
      I1<-max(c(1,(itable[Ips[Iseg],1]-180)))
      Ic<-itable[Ips[Iseg],1]
      I2<-min(c(Nall,(itable[Ips[Iseg],1]+180)))
#     I1<-max(c(1,Ips[Iseg]-180))
#     Ic<-Ips[Iseg]
#     I2<-min(c(N,Ips[Iseg]+180))
      plot(I1:I2,odata[I1:I2,8],type="l",col="black",ylab="",xaxt="n",
           main=paste("Box-Cox transformed dailyP>pthr series\n",
	        "Changepoint ",Iseg," at:",IY0[Ips[Iseg]],sep=""))
      axis(side=1,at=c(I1,Ic,I2),labels=IY0[c(I1,Ic,I2)])
#     tflg<-is.na(p2data[I1:I2,3])==F
#     lines(c(I1:I2)[tflg],p2data[I1:I2,3][tflg],col="red",lwd=1.5)
      lines(c(I1:I2),odata[I1:I2,9],col="red",lwd=1.5)
    }
  }

  yrs<-unique(IY0/10000)
  yrs<-as.integer(seq(yrs[1],yrs[length(yrs)],length=8))
  ymd<-yrs*10000+101
  ats<-rep(NA,length(yrs))
  IY=itable[,2]*10000+itable[,3]*100+itable[,4]
  for(i in 1:length(yrs)){
    it<-match(ymd[i],IY)
    if(!is.na(it)) ats[i]<-itable[it,1]
    else ats[i]<-itable[which.max(IY>ymd[i]),1]
  }
  pdata<-rep(NA,nrow(ori.itable))
  pdata[itable[,1]]<-P
  plot(1:nrow(ori.itable),pdata,type="l",col="black",ylab="prcp(mm)",
       main="Original dailyP>pthr series",xaxt="n")
  axis(side=1,at=ats,labels=yrs)
  lines(1:nrow(ori.itable),odata[,4],col="red")
  if(Ns>0) if(QMout$Mq>1){
    plot(1:nrow(ori.itable),odata[,5],type="l",col="black",ylab="prcp(mm)",
         main="QM adjusted dailyP>pthr series",xaxt="n")
    axis(side=1,at=ats,labels=yrs)
    plot(1:nrow(ori.itable),odata[,6],type="l",col="black",ylab="prcp(mm)",
         main="IBC adjusted dailyP>pthr series",xaxt="n")
    axis(side=1,at=ats,labels=yrs)
  }
  # test plot
  if(Ns>0) if(QMout$Mq>1){
    par(mar=c(4,5,3,2)+.1,cex=.8,mfrow=c(1,1))
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
	       ylab="QM Adjustment",
	       main=paste("distribution of QM adjustments with Mq=",QMout$Mq))
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
        text(.05,ymax-np*(ymax-ymin)/(Ns*3),paste("Seg.",i))
        lines(c(.15,.20),rep(ymax-np*(ymax-ymin)/(Ns*3),2),lwd=2,col=col)
      }
      else np<-np+1
    }
  }
  par(op)
  dev.off()

  otmp<-LSmultiple(Y1,Ti,Ips)
  resi<-otmp$resi
  otmpW<-LSmultiple(W,Ti,Ips)
  resiW<-otmpW$resi
  otmpWL<-LSmultiple(WL,Ti,Ips)
  resiWL<-otmpWL$resi
  otmpWU<-LSmultiple(WU,Ti,Ips)
  resiWU<-otmpWU$resi

  ofileIout<-paste(output,"_fCs.txt",sep="")
  ofileMout<-paste(output,"_mCs.txt",sep="")
  file.create(ofileIout)

  if(Ns==0) {
    cat(paste(Ns,"changepoints in Series","\n"),file=ofileIout)
    cat("PMF finds no Type-1 changepoints in the series!\n")
  }
  else{
    cat(paste(Ns,"changepoints in Series","\n"),
        file=ofileIout)
    for(i in 1:Ns){
      I1<- if(i==1) 1 else Ips[i-1]+1
      Ic<-Ips[i]
      I3<-Ips[i+1]
      Id<-Ids[i]
      Nseg<-I3-I1+1
      PFx95<-getPFx95(cor,Nseg)
      PFx95l<-getPFx95(corL,Nseg)
      PFx95h<-getPFx95(corU,Nseg)
      SSEf.Iseg<-sum(resi[I1:I3]^2)
      Ips1<-Ips[-i]
      otmp1<-LSmultiple(Y1,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      Pk1<-Pk.PMFT(Nseg)
      PFx<-Fx*Pk1[Ic-I1+1]

      SSEf.Iseg<-sum(resiW[I1:I3]^2)
      otmp1<-LSmultiple(W,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(Fx>0) prob<-pf(Fx,1,Nseg-3)
      else{
        Fx<-0
        prob<-0
        PFx<-0
      }

      SSEf.Iseg<-sum(resiWL[I1:I3]^2)
      otmp1<-LSmultiple(WL,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      FxL<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(FxL>0) probL<-pf(Fx,1,Nseg-3)
      else probL<-0

      SSEf.Iseg<-sum(resiWU[I1:I3]^2)
      otmp1<-LSmultiple(WU,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      FxU<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(FxU>0) probU<-pf(Fx,1,Nseg-3)
      else probU<-0
      if(Id==0) { # type-0 changepoints
        if(probU<plev) Idc<-"No  "
	else if(probL<plev&probU>=plev) Idc<-"?   "
	else if(probL>=plev) Idc<-"YifD"
	if(PFx>=PFx95h) Idc<-"Yes "
      }
      else if(Id==1){ # type-1 changepoints
        if(PFx<PFx95l) Idc<-"No  "
        else if(PFx>=PFx95l&PFx<PFx95h) Idc<-"?   "
        else if(PFx>=PFx95h) Idc<-"Yes "
      }

      cat(paste(sprintf("%1.0f",Id)," ",
                sprintf("%-4.4s",Idc),
                sprintf("%10.0f", IY0[Ic])," (",
	        sprintf("%10.4f",probL),"-",
	        sprintf("%10.4f",probU),")",
		sprintf("%6.3f",plev),
	        sprintf("%10.4f",PFx)," (",
	        sprintf("%10.4f",PFx95l),"-",
                sprintf("%10.4f",PFx95h),")\n",sep=""),
		file=ofileIout,
	        append=TRUE)
      cat(paste("PMF : c=", sprintf("%4.0f",Ic), 
      		"; (Time ", sprintf("%10.0f",IY0[Ic]), 
		"); Type= 1; p=",sprintf("%10.4f",prob),"(",
		sprintf("%10.4f",probL),"-",
		sprintf("%10.4f",probU),")",
		"; PFmax=", sprintf("%10.4f",PFx), 
		"; CV95=", sprintf("%10.4f",PFx95), 
		"(", sprintf("%10.4f",PFx95l), 
		"-", sprintf("%10.4f",PFx95h),
		"); Nseg=", sprintf("%4.0f",Nseg),"\n",sep=""), 
		file=ofileSout, append=T)
    }
  }
  if(GUI) return(0) else {
    file.copy(from=ofileIout,to=ofileMout,overwrite=TRUE)
    cat("StepSize.dlyprcp finished successfully...\n")
  }
}
