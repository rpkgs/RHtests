FindU<-function(InSeries,output,MissingValueCode,GUI=FALSE,p.lev=0.95,
                Iadj=10000,Mq=10,Ny4a=0){
  ErrorMSG<-NA
  assign("ErrorMSG",ErrorMSG,envir=.GlobalEnv)
  Nmin<-10
  if(Ny4a>0&Ny4a<=5) Ny4a<-5
  if(!p.lev%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
      ErrorMSG<<-paste("FindU: input p.lev",p.lev,"error\n",
                 get("ErrorMSG",env=.GlobalEnv),"\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
  }
  plev<-p.lev
  pkth<-match(p.lev,c(0.75,0.8,0.9,0.95,0.99,0.9999))
  assign("Nmin",Nmin,envir=.GlobalEnv)
  itmp<-Read(InSeries,MissingValueCode)
  if(itmp==(-1)){
    ErrorMSG<<-paste("FindU: Error in read data from",InSeries,"\n",
               get("ErrorMSG",env=.GlobalEnv),"\n")
    if(!GUI) cat(ErrorMSG)
    return(-1)
  }

  ofileAout<-paste(output,"_U.dat",sep="")
  ofilePdf<-paste(output,"_U.pdf",sep="")
  ofileSout<-paste(output,"_Ustat.txt",sep="")
  file.create(ofileAout)
  file.create(ofilePdf)
  file.create(ofileSout)


  N<-length(Y0); Nadj<-Ny4a*Nt
  cat(paste("The nominal level of confidence (1-alpha)=",plev,"\n"),file=ofileSout)
  cat(paste("Input data filename:", InSeries, "N=",N, "\n"),file=ofileSout,append=T)
  readPFtable(N,pkth)
  Pk0<-Pk.PMFT(N)
  oout<-rmCycle(itable)
  Y1<-oout$Base
  EB<-oout$EB
  assign("EB",EB,envir=.GlobalEnv)
  if(length(EB)!=length(Icy)) {
    ErrorMSG<<-paste("Annual cycle length (from non-missing data) differ from original dataset",
                     "\n",get("ErrorMSG",env=.GlobalEnv),"\n")
    if(GUI==F) print(ErrorMSG)
    return(-1)
  }

  Ip0<-N
  otmp<-LSmultiRedCycle(Y1,Ti,Ip0,1)
  beta0<-otmp$trend
  meanhat0<-otmp$meanhat
  Ehat0<-mean(meanhat0)
  cat(file=ofileSout,paste(" Ignore changepoints -> trend0 =",
      round(beta0,6),"(",round(otmp$betaL,6),",",round(otmp$betaU,6),")",
      "(p=",round(otmp$p.tr,4),"); cor=",
      round(otmp$cor,4),"(", round(otmp$corl,4),",",
      round(otmp$corh,4),")\n"),append=T)

  oout<-PMFT(Y1,Ti,Pk0)
  I0<-0
  I2<-oout$KPx
  I4<-N
  oout1<-PMFxKxI0I2(Y1,Ti,I0,I2)
  I1<-oout1$Ic
  oout2<-PMFxKxI0I2(Y1,Ti,I2,I4)
  I3<-oout2$Ic
  oout3<-PMFxKxI0I2(Y1,Ti,I1,I3)
  I2<-oout3$Ic

  Ns<-1
  Ips<-c(I1,N)
  if(I1>0){
    otmp<-LSmultiple(Y1,Ti,Ips)
    resi<-otmp$resi
    fitted<-otmp$fitted
    otmp<-Rphi(resi,Ips,Ns)
    cor1<-otmp$cor
    corL1<-otmp$corl
    corU1<-otmp$corh
    W<-otmp$W+fitted
    otmp<-PMFxKc(Y1,Ti,I0,I4,I1)
    PFx1<-otmp$PFc
    otmp<-PMFxKc(W,Ti,I0,I4,I1)
    prob1<-otmp$prob
  }
  else{
    prob1<-0
    PFx1<-0
    cor1<-0
    corL1<-0
  }

  Ips<-c(I2,N)
  if(I2>0){
    otmp<-LSmultiple(Y1,Ti,Ips)
    resi<-otmp$resi
    fitted<-otmp$fitted
    otmp<-Rphi(resi,Ips,Ns)
    cor2<-otmp$cor
    corL2<-otmp$corl
    corU2<-otmp$corh
    W<-otmp$W+fitted
    otmp<-PMFxKc(Y1,Ti,I0,I4,I2)
    PFx2<-otmp$PFc
    otmp<-PMFxKc(W,Ti,I0,I4,I2)
    prob2<-otmp$prob
  }
  else{
    prob2<-0
    PFx2<-0
    cor2<-0
    corL2<-0
  }

  Ips<-c(I3,N)
  if(I3>0){
    otmp<-LSmultiple(Y1,Ti,Ips)
    resi<-otmp$resi
    fitted<-otmp$fitted
    otmp<-Rphi(resi,Ips,Ns)
    cor3<-otmp$cor
    corL3<-otmp$corl
    corU3<-otmp$corh
    W<-otmp$W+fitted
    otmp<-PMFxKc(Y1,Ti,I0,I4,I3)
    PFx3<-otmp$PFc
    otmp<-PMFxKc(W,Ti,I0,I4,I3)
    prob3<-otmp$prob
  }
  else{
    prob3<-0
    PFx3<-0
    cor3<-0
    corL3<-0
  }
 
  ofileIout<-paste(output,"_1Cs.txt",sep="")
  ofileMout<-paste(output,"_mCs.txt",sep="")
  file.create(ofileIout)

  tmp<-sort(c(PFx1,PFx2,PFx3),decreasing=T,index.return=T)
  PFx.mx<-tmp$x[1]
  prob.mx<-c(prob1,prob2,prob3)[tmp$ix[1]]
  Imx<-c(I1,I2,I3)[tmp$ix[1]]
  cor.mx<-c(corL1,corL2,corL3)[tmp$ix[1]]
  PFx95L<-getPFx95(cor.mx,N)
  if(PFx.mx<PFx95L){
    Ns<-0
    Ips<-N
    cat(paste(Ns,"changepoints in Series", InSeries,
              paste("sample:(",sprintf("%1.0f",1)," ",sprintf("%-4.4s","YifD"),
	       sprintf("%10.0f",19000101),")",sep=""),"\n"),file=ofileIout)
    cat("PMF finds no Type-1 changepoints in the series!\n")
#   return(0)
  }

  else{
  Ns<-1
  Ips<-c(Imx,N)

# if(Debug) cat(file=flog,c("First Ips:",Ips,"\n"))
  tt<-TRUE
  Niter<-0
  while(tt){  # condition on is there any more Bps to insert in Ips?
    Niter<-Niter+1
    tt<-FALSE
    Ips0<-NULL
    for(i in 1:(Ns+1)){ # search for new break points
      I0<- if(i==1) 0 else Ips[i-1]
      I2<-Ips[i]
      otmp<-PMFxKxI0I2(Y1,Ti,I0,I2)
      if(otmp$prob>0) Ips0<-sort(c(Ips0,otmp$Ic))
    }
#   if(Debug) {
#     cat(file=flog,paste("Niter",Niter,"new Ips:",length(Ips0),"\n"),append=T)
#     cat(file=flog,c(Ips0,"\n"),append=T)
#     cat(file=flog,paste("Niter",Niter,"old Ns",Ns," "),append=T)
#   }
 # finished search for possible new changepoints, start estimate the p-value
 # of the Ips0 series and find the most significant changepoint Iseg.mx
    tt1<-TRUE
    while(tt1){
      PFx.mx<-(-9999)
      Iseg.mx<-0
      PFx95L.mx<-0
      if(length(Ips0)==0) tt1<-FALSE
      else{
        for(i in 1:length(Ips0)){
	  Ips1<-sort(c(Ips,Ips0[i]))
	  ith<-match(Ips0[i],Ips1)
	  otmp<-PMFxIseg(Y1,Ti,Ips1,ith)
	  if(otmp$PFx<otmp$PFx95L) Ips0[i]<-0
	  else
	    if(otmp$PFx>PFx.mx){
	      PFx.mx<-otmp$PFx
	      Iseg.mx<-Ips0[i]
	      PFx95L.mx<-otmp$PFx95L
	    }
	}
	if(PFx.mx>=PFx95L.mx){
	  Ips<-sort(c(Ips,Iseg.mx))
	  Ns<-Ns+1
	  Ips0<-Ips0[Ips0!=Iseg.mx]
	  tt<-TRUE
	}
	else tt1<-FALSE
	Ips0<-Ips0[Ips0!=0]
      }
    }
#   if(Debug) cat(file=flog,paste("new Ns:", Ns,"\n"),append=T)
  }  
# finish finding any possible new changepoints
# start to delete the least significant changepoint if it is insignificant

  tt<-TRUE
  while(tt){
    tt<-FALSE
    Iseg.mn<-0
    PFx.mn<-99999
    PFx95L.mn<-99999
    for(i in 1:Ns){
      otmp<-PMFxIseg(Y1,Ti,Ips,i)
      if(otmp$PFx<PFx.mn){
        Iseg.mn<-i
	PFx.mn<-otmp$PFx
	PFx95L.mn<-otmp$PFx95L
      }
    }
    if(Iseg.mn>0&PFx.mn<PFx95L.mn){
      Ips<-Ips[-Iseg.mn]
      Ns<-Ns-1
      if(Ns>0) tt<-TRUE
    }
  }

  }
# end of detection
  if(Ns>0) {
    Nsegs<-Ips-c(0,Ips[1:Ns])
    Iseg.longest<-sort(Nsegs,index=T,decreasing=T)$ix[1]
  }
  else Iseg.longest<-0

  if(Iadj>(Ns+1)|Iseg.longest==0) Iseg.adj<-Ns+1 
  else if(Iadj==0)Iseg.adj<-Iseg.longest
  else Iseg.adj<-Iadj

  oout<-LSmultiRedCycle(Y1,Ti,Ips,Iseg.adj)
  Y1<-oout$Y0
  cor<-oout$cor
  corl<-oout$corl
  corh<-oout$corh
  df<-(N-2-Nt-Ns)
  pcor<-pt(abs(cor)*sqrt(df/(1-cor^2)),df)
  W<-oout$W
  WL<-oout$WL
  WU<-oout$WU
  EB1<-oout$EB
  itmp1<-cbind(EB1,Icy)
  itmp2<-cbind(1:N,Imd)
  colnames(itmp2)<-c("idx","Icy")
  itmp<-merge(itmp1,itmp2,by="Icy")
  EBfull<-itmp[order(itmp[,"idx"]),"EB1"]
  EEB<-mean(EB1,na.rm=T)

  if(Ns>0){
    Rb<-Y1-oout$trend*Ti+EBfull
    QMout<-QMadjGaussian(Rb,Ips,Mq,Iseg.adj,Nadj)
    B<-QMout$PA
    cat(paste("Nseg_shortest =",QMout$Nseg.mn,"; Mq = ",QMout$Mq,"; Ny4a = ",Ny4a,"\n"),
        file=ofileSout,append=T)
    cat(paste("\n Adjust to segment", Iseg.adj,": from",
        if(Iseg.adj==1) 1 else Ips[Iseg.adj-1]+1,
        "to",Ips[Iseg.adj],"\n"),file=ofileSout,append=T)
#   cat("#Fcat, DP (CDF and Differnces in category mean)\n",file=ofileSout,
#       append=T)
    if(Mq>1){
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
        cat(paste("Seg. ",i,": mean of QM-adjustments =",round(mean(QMout$W[I1:I2]),4),
            "\n",sep=""),file=ofileSout,append=T)
      }
    }
  }
  else B<-Y1-oout$trend*Ti+EBfull

# itmp1<-cbind(EB1,Icy)
# itmp2<-cbind(1:N,Imd)
# colnames(itmp2)<-c("idx","Icy")
# itmp<-merge(itmp1,itmp2,by="Icy")
# EBfull<-itmp[order(itmp[,"idx"]),"EB1"]
# EEB<-mean(EB1,na.rm=T)
  adj<-oout$Y0+EBfull
# B<-B+EBfull+oout$trend*Ti
  B<-B+oout$trend*Ti
  cat("Common trend TPR fit to the de-seasonalized Base series:\n",
      file=ofileSout,append=T)
  cat(paste("#steps= ",Ns,"; trend=",round(oout$trend,6),"(",
            round(oout$betaL,6),",",round(oout$betaU,6),") (p=",
            round(oout$p.tr,4),"); cor=",
            round(cor,4),"(",round(corl,4),",",round(corh,4),")",
	    round(pcor,4),"\n"),
	    file=ofileSout,append=T)
  if(Ns>0) for(i in 1:(Ns+1)){
    I1<-if(i==1) 1 else Ips[i-1]+1
    I2<-Ips[i]
    Delta<-oout$mu[Iseg.adj]-oout$mu[i]
    adj[I1:I2]<-adj[I1:I2]+Delta
    stepsize<-oout$mu[i+1]-oout$mu[i]
    cat(paste(Ips[i],IY0[Ips[i]],"stepsize=",round(stepsize,4),"\n"),
        file=ofileSout,append=T)
  }

  oR<-Y1-oout$meanhat
  oR[2:N]<-oR[2:N]-oR[1:(N-1)]*cor
  Ehat<-mean(oout$meanhat)
  meanhat0<-meanhat0-Ehat0+Ehat

  pdf(file=ofilePdf,onefile=T,paper='letter')
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  par(mfrow=c(2,1))
  par(mar=c(3,4,3,2)+.1,cex.main=.8,cex.lab=.8,cex.axis=.8,cex=.8)
  uyrs<-unique(floor(ori.itable[,1]/10))*10
  labels<-NULL
  ats<-NULL
  for(i in 1:length(uyrs)){
    if(!is.na(match(uyrs[i],ori.itable[,1]))){
      labels<-c(labels,uyrs[i])
      ats<-c(ats,match(uyrs[i],ori.itable[,1]))
    }
  }
  
  pdata<-rep(NA,dim(ori.itable)[1])
  pdata[ooflg]<-oout$Y0
  plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(oout$Y0,oout$meanhat),max(oout$Y0,oout$meanhat)),
       xaxt="n",col="black",lwd=0.5,
       main="Base anomaly series and regression fit")
  axis(side=1,at=ats,labels=labels)
  pdata[ooflg]<-oout$meanhat
  lines(1:dim(ori.itable)[1],pdata,col="red")

  pdata[ooflg]<-oout$Y0+EBfull
  plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(oout$Y0+EBfull,oout$meanhat+EEB),
              max(oout$Y0+EBfull,oout$meanhat+EEB)),
       xaxt="n",col="black",lwd=0.5,
       main="Base series and regression fit")
  axis(side=1,at=ats,labels=labels)
  pdata[ooflg]<-oout$meanhat+EEB
  lines(1:dim(ori.itable)[1],pdata,col="red")

  pdata[ooflg]<-adj
  plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(c(adj,B)),max(c(adj,B))),
       xaxt="n",col="black",lwd=0.5,
       main="Mean-adjusted base series")
  axis(side=1,at=ats,labels=labels)

  if(Mq>1){
    pdata[ooflg]<-B
    plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
         ylim=c(min(c(adj,B)),max(c(adj,B))),
         xaxt="n",col="black",lwd=0.5,
         main="QM-adjusted base series")
    axis(side=1,at=ats,labels=labels)
  }

  # test plot
  if(Ns>0&Mq>1){
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
          title(cex.main=0.9,main=paste("distribution of QM adjustments with Mq=",QMout$Mq),line=.5)
	  icol<-2*np
	  for(j in 1:QMout$Mq){
	    lines(c(osmean[(j+1),icol]-Fd,osmean[(j+1),icol]+Fd),
	          c(rep(osmean[(j+1),(icol+1)],2)),col=col,lty=2,lwd=.5)
	    if(j>=1&j<QMout$Mq) lines(rep(osmean[(j+1),icol]+Fd,2),
	          c(osmean[(j+1),(icol+1)],osmean[(j+2),(icol+1)]),col=col,lty=2,lwd=0.5)
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
	          c(osmean[(j+1),(icol+1)],osmean[(j+2),(icol+1)]),col=col,lty=2,lwd=0.5)
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

  odata<-matrix(NA,dim(ori.itable)[1],10)
# odata[ooflg,1]<-Ti
  odata[,1]<-c(1:dim(ori.itable)[1])
  odata[,2]<-ori.itable[,1]*10000+ori.itable[,2]*100+ori.itable[,3]
# odata[ooflg,3]<-round(oout$Y0+EBfull,4)
  odata[,3]<-ori.itable[,4]
  odata[ooflg,4]<-round(oout$meanhat+EEB,4)
  odata[ooflg,5]<-round(adj,4)
  odata[ooflg,6]<-round(oout$Y0,4)
  odata[ooflg,7]<-round(oout$meanhat,4)
  odata[ooflg,8]<-round(oout$meanhat+EBfull,4)
  if(Ns>0) if(QMout$Mq>1) odata[ooflg,9]<-round(B,4)
  odata[ooflg,10]<-round(meanhat0,4)

  Imd1<-ori.itable[,2]*100+ori.itable[,3]
  if(sum(is.na(ori.itable[,4])==F&Imd1==229)>0){
    if(Ns>0){
      tdata<-ori.itable[is.na(ori.itable[,4])==F,]
      IY1<-tdata[,1]*10000+tdata[,2]*100+tdata[,3]
      Ips.ymd<-IY0[Ips]
      Ips.1<-rep(NA,Ns+1)
      for(i in 1:Ns) Ips.1[i]<-c(1:length(IY1))[IY1==Ips.ymd[i]]
      Ips.1[Ns+1]<-length(IY1)
#     Ips.1<-c(1:length(IY1))[Ips.ymd==IY1]
      Imd2<-tdata[,2]*100+tdata[,3]
      Ids.leap<-c(1:length(Imd2))[Imd2==229]
      Nl<-length(Ids.leap)
      Rb<-Y1-oout$trend*Ti+EBfull
      Rb1<-tdata[,4]; Rb1[-Ids.leap]<-Rb
      Ti1<-rep(NA,length(IY1)); Ti1[-Ids.leap]<-Ti
      for(i in 1:length(Ids.leap)) {
        Rb1[Ids.leap[i]]<-tdata[Ids.leap[i],4]+Rb1[Ids.leap[i]-1]-tdata[Ids.leap[i]-1,4]
        Ti1[Ids.leap[i]]<-Ti1[Ids.leap[i]-1]
      }
      if(QMout$Mq>1){
        B1<-QMadjGaussian(Rb1,Ips.1,Mq,Iseg.adj,Nadj)$PA
        B1<-B1+oout$trend*Ti1
        B1.leap<-B1[Ids.leap]
        odata[is.na(odata[,3])==F&Imd1==229,9]<-round(B1.leap,4)
      }
    }
    else
      odata[Imd1==229,9]<-odata[Imd1==229,3]
    Ids.leapo<-c(1:dim(ori.itable)[1])[is.na(ori.itable[,4])==F&Imd1==229]
    for(jth in 1:length(Ids.leapo)){
      kth<-Ids.leapo[jth]
      if(Ns>0){
        k1th<-if(odata[kth-1,2]%in%IY0[Ips]) (kth+1) else (kth-1)
      }
      else k1th<-kth-1
      for(pth in c(4,7,8,10)) odata[kth,pth]<-odata[k1th,pth]
      for(pth in c(5,6)){delta1<-odata[k1th,3]-odata[k1th,pth]; odata[kth,pth]<-odata[kth,3]-delta1}
    }
  }
    
  write.table(file=ofileAout,odata,na=MissingValueCode,
	col.names=F,row.names=F)

  otmp<-LSmultiple(Y1,Ti,Ips)
  resi<-otmp$resi
  otmpW<-LSmultiple(W,Ti,Ips)
  resiW<-otmpW$resi
  otmpWL<-LSmultiple(WL,Ti,Ips)
  resiWL<-otmpWL$resi
  otmpWU<-LSmultiple(WU,Ti,Ips)
  resiWU<-otmpWU$resi
  if(Ns==0) {
    cat(paste(Ns,"changepoints in Series", InSeries,
              paste("sample:(",sprintf("%1.0f",1)," ",sprintf("%-4.4s","YifD"),
	       sprintf("%10.0f",19000101),")",sep=""),"\n"),file=ofileIout)
    cat("PMF finds no Type-1 changepoints in the series!\n")
#   return(0)
  }
  else{
    cat(paste(Ns,"changepoints in Series", InSeries,"\n"),
        file=ofileIout)
    for(i in 1:Ns){
      I1<- if(i==1) 1 else Ips[i-1]+1
      Ic<-Ips[i]
      I3<-Ips[i+1]
      Nseg<-I3-I1+1
      PFx95<-getPFx95(cor,Nseg)
      PFx95l<-getPFx95(corl,Nseg)
      PFx95h<-getPFx95(corh,Nseg)
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

#     else if(Id==1) { # type-1 changepoints
        if(PFx<PFx95l) Idc<-"No  "
        else if(PFx>=PFx95l&PFx<PFx95h) Idc<-"?   "
        else if(PFx>=PFx95h) Idc<-"Yes "
#     }
      cat(paste(sprintf("%1.0f",1)," ",
                sprintf("%-4.4s",Idc),
                sprintf("%10.0f", IY0[Ic])," (",
	        sprintf("%6.4f",probL),"-",
	        sprintf("%6.4f",probU),")",
		sprintf("%6.3f",plev),
	        sprintf("%10.4f",PFx)," (",
	        sprintf("%10.4f",PFx95l),"-",
                sprintf("%10.4f",PFx95h),")\n",sep=""),
		file=ofileIout,
	        append=TRUE)
      cat(paste("PMF : c=", sprintf("%4.0f",Ic), 
      		"; (Time ", sprintf("%10.0f",IY0[Ic]), 
		"); Type= 1; p=",sprintf("%10.4f",prob),"(",
		sprintf("%6.4f",probL),"-",
		sprintf("%6.4f",probU),")",
		"; PFmax=", sprintf("%10.4f",PFx), 
		"; CV95=", sprintf("%10.4f",PFx95), 
		"(", sprintf("%10.4f",PFx95l), 
		"-", sprintf("%10.4f",PFx95h),
		"); Nseg=", sprintf("%4.0f",Nseg),"\n",sep=""), 
		file=ofileSout, append=T)
    }
  }
  if(GUI)
    return(0)
  else{
    file.copy(from=ofileIout,to=ofileMout,overwrite=TRUE)
    cat("FindU finished successfully...\n")
  }
}

FindUD<-function(InSeries,InCs,output,MissingValueCode,GUI=FALSE,p.lev=0.95,
                 Iadj=10000,Mq=10,Ny4a=0){
  Debug<-TRUE
  ErrorMSG<-NA
  assign("ErrorMSG",ErrorMSG,envir=.GlobalEnv)
  flog<-paste(output,".log",sep="")
  Nmin<-10
  if(Ny4a>0&Ny4a<=5) Ny4a<-5
  if(!p.lev%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
      ErrorMSG<<-paste("FindU: input p.lev",p.lev,"error\n",
                 get("ErrorMSG",env=.GlobalEnv),"\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
  }
  plev<-p.lev
  pkth<-match(p.lev,c(0.75,0.8,0.9,0.95,0.99,0.9999))
  assign("Nmin",Nmin,envir=.GlobalEnv)
  itmp<-Read(InSeries,MissingValueCode)
  if(itmp<0){
    ErrorMSG<<-paste("FindUD: Error in read data from",InSeries,"\n",
               get("ErrorMSG",env=.GlobalEnv),"\n")
    if(!GUI) cat(ErrorMSG)
    return(-1)
  }
  N<-length(Y0); Nadj<-Ny4a*Nt
  readPFtable(N,pkth)
  itmp<-readLines(InCs)
  Pk0<-Pk.PMFT(N)
  oout<-rmCycle(itable)
  Y1<-oout$Base
  EB<-oout$EB
  assign("EB",EB,envir=.GlobalEnv)
  if(length(EB)!=length(Icy)) {
    ErrorMSG<<-paste("Annual cycle length (from non-missing data) differ from original dataset",
                     "\n",get("ErrorMSG",env=.GlobalEnv),"\n")
    if(!GUI) cat(ErrorMSG)
    return(-1)
  }

  ofileIout<-paste(output,"_pCs.txt",sep="")
  ofileMout<-paste(output,"_mCs.txt",sep="")
  file.create(ofileIout)

  ofileAout<-paste(output,"_UD.dat",sep="")
  ofilePdf<-paste(output,"_UD.pdf",sep="")
  ofileSout<-paste(output,"_UDstat.txt",sep="")
  file.create(ofileAout)
  file.create(ofilePdf)
  file.create(ofileSout)
  cat(paste("The nominal level of confidence (1-alpha)=",plev,"\n"),file=ofileSout)
  cat(paste("Input data filename:", InSeries,"N=",N,"\n"),file=ofileSout,append=T)

  Ip0<-N
  oout<-LSmultiRedCycle(Y1,Ti,Ip0,1)
  beta0<-oout$trend
  betaL0<-oout$betaL
  betaU0<-oout$betaU
  meanhat0<-oout$meanhat
  Ehat0<-mean(meanhat0)
  p.tr0<-oout$p.tr
  corD<-oout$cor
  corDL<-oout$corl
  corDU<-oout$corh

  if(length(itmp)<2){ # no input changepoints
    oout<-PMFT(Y1,Ti,Pk0)
    I0<-0
    I2<-oout$KPx
    I4<-N
    oout1<-PMFxKxI0I2(Y1,Ti,I0,I2)
    I1<-oout1$Ic
    oout2<-PMFxKxI0I2(Y1,Ti,I2,I4)
    I3<-oout2$Ic
    oout3<-PMFxKxI0I2(Y1,Ti,I1,I3)
    I2<-oout3$Ic

    Ns<-1
    Ips<-c(I1,N)
    if(I1>0){
      otmp<-LSmultiple(Y1,Ti,Ips)
      resi<-otmp$resi
      fitted<-otmp$fitted
      otmp<-Rphi(resi,Ips,Ns)
      cor1<-otmp$cor
      corL1<-otmp$corl
      corU1<-otmp$corh
      W<-otmp$W+fitted
      otmp<-PMFxKc(Y1,Ti,I0,I4,I1)
      PFx1<-otmp$PFc
      otmp<-PMFxKc(W,Ti,I0,I4,I1)
      prob1<-otmp$prob
    }
    else{
      prob1<-0
      PFx1<-0
    }

    Ips<-c(I2,N)
    if(I2>0){
      otmp<-LSmultiple(Y1,Ti,Ips)
      resi<-otmp$resi
      fitted<-otmp$fitted
      otmp<-Rphi(resi,Ips,Ns)
      cor2<-otmp$cor
      corL2<-otmp$corl
      corU2<-otmp$corh
      W<-otmp$W+fitted
      otmp<-PMFxKc(Y1,Ti,I0,I4,I2)
      PFx2<-otmp$PFc
      otmp<-PMFxKc(W,Ti,I0,I4,I2)
      prob2<-otmp$prob
    }
    else{
      prob2<-0
      PFx2<-0
    }

    Ips<-c(I3,N)
    if(I3>0){
      otmp<-LSmultiple(Y1,Ti,Ips)
      resi<-otmp$resi
      fitted<-otmp$fitted
      otmp<-Rphi(resi,Ips,Ns)
      cor3<-otmp$cor
      corL3<-otmp$corl
      corU3<-otmp$corh
      W<-otmp$W+fitted
      otmp<-PMFxKc(Y1,Ti,I0,I4,I3)
      PFx3<-otmp$PFc
      otmp<-PMFxKc(W,Ti,I0,I4,I3)
      prob3<-otmp$prob
    }
    else{
      prob3<-0
      PFx3<-0
    }
  
    tmp<-sort(c(PFx1,PFx2,PFx3),decreasing=T,index.return=T)
    PFx.mx<-tmp$x[1]
    prob.mx<-c(prob1,prob2,prob3)[tmp$ix[1]]
    Imx<-c(I1,I2,I3)[tmp$ix[1]]
    if(prob.mx<plev){
#     cat("PMF finds the series to be homogeneous!\n",file=ofileIout)
      cat(paste(0,"changepoints in Series", InSeries,"\n"),file=ofileIout)
      cat("PMF finds the series to be homogeneous!\n")
#     return()
      Ns<-0
      Ips<-N
      Ids<-c(0)
    }
    else{
      Ns<-1
      Ips<-c(Imx,N)
      Ids<-c(0,1)
    }
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
      ErrorMSG<<-paste("FindUD: Ips read in from ",InCs,"error:")
	for(i in 1:Ns)
        ErrorMSG<<-paste(get("ErrorMSG",env=.GlobalEnv),Ips[i])
      ErrorMSG<<-paste(get("ErrorMSG",env=.GlobalEnv),"\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
    }
  }

  if(Ns>0){
  Ips.i<-Ips
  Niter<-0
# start search for all possible changepoints
  tt<-TRUE
  while(tt){
    Niter<-Niter+1
    tt<-FALSE
    Ips0<-NULL
    for(i in 1:(Ns+1)){
      I0<- if(i==1) 0 else Ips[i-1]
      I2<-Ips[i]
      otmp<-PMFxKxI0I2(Y1,Ti,I0,I2)
      if(otmp$prob>0) Ips0<-sort(c(Ips0,otmp$Ic))
    }
# estimate p-value of each changepoint in series Ips0, find the most significant
# changepoint Iseg.mx
    tt1<-TRUE
    while(tt1){
      if(length(Ips0)==0) tt1<-FALSE
      else{
        Iseg.mx<-0
	prob.mx<-(-1)
	probL.mx<-(-1)
	PFx.mx<-(-1)
        for(i in 1:length(Ips0)){
	  Ips1<-sort(c(Ips,Ips0[i]))
	  ith<-match(Ips0[i],Ips1)
	  otmp<-PMFxIseg(Y1,Ti,Ips1,ith)
	  probL<-min(c(otmp$probL,otmp$probU,otmp$prob))
	  probU<-max(c(otmp$probL,otmp$probU,otmp$prob))
	  PFx<-otmp$PFx
	  if(probU<plev) Ips0[i]<-0
	  else
	    if(PFx>PFx.mx){
	      prob.mx<-otmp$prob
	      probL.mx<-probL
	      Iseg.mx<-Ips0[i]
	      PFx.mx<-PFx
	    }
	}
	if(probL.mx>=plev){
	  Ips<-sort(c(Ips,Iseg.mx))
	  Ns<-Ns+1
	  Ips0<-Ips0[Ips0!=Iseg.mx]
	  tt<-TRUE
	}
	else tt1<-FALSE
	Ips0<-Ips0[Ips0!=0]
      }
    }
  }
  Ids0<-rep(NA,length(Ips))
  for(i in 1:length(Ips)){
    if(Ips[i]%in%Ips.i) Ids0[i]<-Ids[Ips.i==Ips[i]]
    else Ids0[i]<-0
  }
  Ids<-Ids0
# Ids<-as.integer(Ips%in%Ips.i)
  tt<-TRUE
  while(tt){
    tt<-FALSE
    probL.mn<-9999
    Iseg.mn<-0
    for(i in 1:Ns){
      if(Ids[i]==0){ # check those un-documented
        Ips0<-Ips[-i]
        otmp<-PMFxIseg(Y1,Ti,Ips,i)
        probL<-min(otmp$probL,otmp$probU)
        if(probL<probL.mn){
          Iseg.mn<-i
          probL.mn<-probL
        }
      } # end if documented
    }
    if(Iseg.mn>0&probL.mn<plev){
      Ips<-Ips[-Iseg.mn]
      Ids<-Ids[-Iseg.mn]
      Ns<-Ns-1
      if(Ns>0) tt<-TRUE
    }
  }
  }
# all changepoints are significant, calculate stats and output

  if(Ns>0) {
    Nsegs<-Ips-c(0,Ips[1:Ns])
    Iseg.longest<-sort(Nsegs,index=T,decreasing=T)$ix[1]
  }
  else Iseg.longest<-0

  if(Iadj>(Ns+1)|Iseg.longest==0) Iseg.adj<-Ns+1 
  else if(Iadj==0)Iseg.adj<-Iseg.longest
  else Iseg.adj<-Iadj

  otmp<-LSmultiRedCycle(Y1,Ti,Ips,Iseg.adj)
  Y1<-otmp$Y0
  cor<-otmp$cor
  corl<-otmp$corl
  corh<-otmp$corh
  df<-(N-2-Nt-Ns)
  pcor<-pt(abs(cor)*sqrt(df/(1-cor^2)),df)
  Rf<-otmp$resi
  W<-otmp$W
  WL<-otmp$WL
  WU<-otmp$WU
  EB1<-otmp$EB

  itmp1<-cbind(EB1,Icy)
  itmp2<-cbind(1:N,Imd)
  colnames(itmp2)<-c("idx","Icy")
  itmp<-merge(itmp1,itmp2,by="Icy")
  EBfull<-itmp[order(itmp[,"idx"]),"EB1"]
  EEB<-mean(EB1,na.rm=T)

  if(Ns>0){
    Rb<-Y1-otmp$trend*Ti+EBfull
    QMout<-QMadjGaussian(Rb,Ips,Mq,Iseg.adj,Nadj)
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
        cat(paste("Seg. ",i,": mean of QM-adjustments =",round(mean(QMout$W[I1:I2]),4),
            "\n",sep=""),file=ofileSout,append=T)
      }
    }
  }
  else B<-Y1-otmp$trend*Ti+EBfull

  adj<-otmp$Y0+EBfull
# B<-B+EBfull+otmp$trend*Ti
  B<-B+otmp$trend*Ti

  cat(file=ofileSout,paste(" Ignore changepoints -> trend0 =",
      round(beta0,6),"(",round(betaL0,6),",",round(betaU0,6),") (p=",
      round(p.tr0,4),"); cor=", round(corD,4),"(",round(corDL,4),
      ",",round(corDU,4),")\n"),append=T)
  cat("Common trend TPR fit to the seasonal-cycle-adjusted Base series:\n",
      file=ofileSout,append=T)
  cat(paste("#steps= ",Ns,"; trend=",round(otmp$trend,6),"(",
            round(otmp$betaL,6),",",round(otmp$betaU,6),") (p=",
            round(otmp$p.tr,4),"); cor=",
            round(cor,4),"(",round(corl,4),",",round(corh,4),")",
	    round(pcor,4),"\n"),
	    file=ofileSout,append=T)
  if(Ns>0) for(i in 1:(Ns+1)){
    I1<-if(i==1) 1 else Ips[i-1]+1
    I2<-Ips[i]
    Delta<-otmp$mu[Iseg.adj]-otmp$mu[i]
    adj[I1:I2]<-adj[I1:I2]+Delta
    stepsize<-otmp$mu[i+1]-otmp$mu[i]
    cat(paste(Ips[i],IY0[Ips[i]],"stepsize=",round(stepsize,4),"\n"),
        file=ofileSout,append=T)
  }

  oR<-Y1-otmp$meanhat
  oR[2:N]<-oR[2:N]-oR[1:(N-1)]*cor
  Ehat<-mean(otmp$meanhat)
  meanhat0<-meanhat0-Ehat0+Ehat

  pdf(file=ofilePdf,onefile=T,paper='letter')
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  par(mfrow=c(2,1))
  par(mar=c(3,4,3,2)+.1,cex.main=.8,cex.lab=.8,cex.axis=.8,cex=.8)

  uyrs<-unique(floor(ori.itable[,1]/10))*10
  labels<-NULL
  ats<-NULL
  for(i in 1:length(uyrs)){
    if(!is.na(match(uyrs[i],ori.itable[,1]))){
      labels<-c(labels,uyrs[i])
      ats<-c(ats,match(uyrs[i],ori.itable[,1]))
    }
  }

  pdata<-rep(NA,dim(ori.itable)[1])
  pdata[ooflg]<-otmp$Y0
  plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(otmp$Y0,otmp$meanhat),max(otmp$Y0,otmp$meanhat)),
       xaxt="n",col="black",lwd=0.5,
       main="Base anomaly series and regression fit")
  axis(side=1,at=ats,labels=labels)
  pdata[ooflg]<-otmp$meanhat
  lines(1:dim(ori.itable)[1],pdata,col="red")

  pdata[ooflg]<-otmp$Y0+EBfull
  plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(otmp$Y0+EBfull,otmp$meanhat+EEB),
              max(otmp$Y0+EBfull,otmp$meanhat+EEB)),
       xaxt="n",col="black",lwd=0.5,
       main="Base series and regression fit")
  axis(side=1,at=ats,labels=labels)

  pdata[ooflg]<-otmp$meanhat+EEB
  lines(1:dim(ori.itable)[1],pdata,col="red")
  
  pdata[ooflg]<-adj
  plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(c(adj,B)),max(c(adj,B))),
       xaxt="n",col="black",lwd=0.5,
       main="Mean-adjusted base series")
  axis(side=1,at=ats,labels=labels)

  if(QMout$Mq>1){
    pdata[ooflg]<-B
    plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
         ylim=c(min(c(adj,B)),max(c(adj,B))),
         xaxt="n",col="black",lwd=0.5,
         main="QM-adjusted base series")
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
          title(cex.main=0.9,main=paste("distribution of QM adjustments with Mq=",QMout$Mq),line=0.5)
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

  odata<-matrix(NA,dim(ori.itable)[1],10)
# odata[ooflg,1]<-Ti
  odata[,1]<-c(1:dim(ori.itable)[1])
  odata[,2]<-ori.itable[,1]*10000+ori.itable[,2]*100+ori.itable[,3]
# odata[ooflg,3]<-round(otmp$Y0+EBfull,4)
  odata[,3]<-ori.itable[,4]
  odata[ooflg,4]<-round(otmp$meanhat+EEB,4)
  odata[ooflg,5]<-round(adj,4)
  odata[ooflg,6]<-round(otmp$Y0,4)
  odata[ooflg,7]<-round(otmp$meanhat,4)
  odata[ooflg,8]<-round(otmp$meanhat+EBfull,4)
  if(Ns>0) if(QMout$Mq>1) odata[ooflg,9]<-round(B,4)
# odata[ooflg,9]<-round(oR,4)
# odata[ooflg,10]<-round(Rb,4)
  odata[ooflg,10]<-round(meanhat0,4)

  Imd1<-ori.itable[,2]*100+ori.itable[,3]
  if(sum(is.na(ori.itable[,4])==F&Imd1==229)>0){
    if(Ns>0){
      tdata<-ori.itable[is.na(ori.itable[,4])==F,]
      IY1<-tdata[,1]*10000+tdata[,2]*100+tdata[,3]
      Ips.ymd<-IY0[Ips]
      Ips.1<-rep(NA,Ns+1)
      for(i in 1:Ns) Ips.1[i]<-c(1:length(IY1))[IY1==Ips.ymd[i]]
      Ips.1[Ns+1]<-length(IY1)
#     Ips.1<-c(1:length(IY1))[Ips.ymd==IY1]
      Imd2<-tdata[,2]*100+tdata[,3]
      Ids.leap<-c(1:length(Imd2))[Imd2==229]
      Nl<-length(Ids.leap)
      Rb<-Y1-otmp$trend*Ti+EBfull
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
        odata[is.na(odata[,3])==F&Imd1==229,9]<-round(B1.leap,4)
      }
    }
    else
      odata[Imd1==229,9]<-odata[Imd1==229,3]
    Ids.leapo<-c(1:dim(ori.itable)[1])[is.na(ori.itable[,4])==F&Imd1==229]
    for(jth in 1:length(Ids.leapo)){
      kth<-Ids.leapo[jth]
      if(Ns>0){
        k1th<-if(odata[kth-1,2]%in%IY0[Ips]) (kth+1) else (kth-1)
      }
      else k1th<-kth-1
      for(pth in c(4,7,8,10)) odata[kth,pth]<-odata[k1th,pth]
      for(pth in c(5,6)){delta1<-odata[k1th,3]-odata[k1th,pth]; odata[kth,pth]<-odata[kth,3]-delta1}
    }
  }
    
  write.table(file=ofileAout,odata,na=MissingValueCode,
	col.names=F,row.names=F)
  otmp<-LSmultiple(W,Ti,Ips)
  RW<-otmp$resi
  otmp<-LSmultiple(WL,Ti,Ips)
  RWL<-otmp$resi
  otmp<-LSmultiple(WU,Ti,Ips)
  RWU<-otmp$resi
  if(Ns==0) {
#   cat("PMF finds the series to be homogeneous!\n",file=ofileIout)
    cat(paste(Ns,"changepoints in Series", InSeries,"\n"),file=ofileIout)
    cat("PMF finds the series to be homogeneous!\n")
#   return()
  }
  else{
    cat(paste(Ns,"changepoints in Series", InSeries,"\n"),
        file=ofileIout)

    for(i in 1:Ns){
      I1<- if(i==1) 1 else Ips[i-1]+1
      I3<-Ips[i+1]
      Ic<-Ips[i]
      Id<-Ids[i]
      Nseg<-I3-I1+1
    
      PFx95<-getPFx95(cor,Nseg)
      PFx95L<-getPFx95(corl,Nseg)
      PFx95U<-getPFx95(corh,Nseg)
      SSEf.Iseg<-sum(Rf[I1:I3]^2)
      Ips0<-Ips[-i]
      otmp<-LSmultiple(Y1,Ti,Ips0)
      SSE0.Iseg<-sum(otmp$resi[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      Pk0<-Pk.PMFT(Nseg)
      PFx<-Fx*Pk0[Ic-I1+1]

      otmp<-LSmultiple(W,Ti,Ips0)
      SSE0.Iseg<-sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg<-sum(RW[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(Fx<=0){
        PFx<-0
        Fx<-0
        prob<-0
      }
      else prob<-pf(Fx,1,Nseg-3)

      otmp<-LSmultiple(WL,Ti,Ips0)
      SSE0.Iseg<-sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg<-sum(RWL[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      probL0<-if(Fx<0) 0 else pf(Fx,1,Nseg-3)

      otmp<-LSmultiple(WU,Ti,Ips0)
      SSE0.Iseg<-sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg<-sum(RWU[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      probU0<-if(Fx<0) 0 else pf(Fx,1,Nseg-3)

      probL<-min(probL0,probU0)
      probU<-max(probL0,probU0)

      if(Id==0) { # type-0 changepoints
        if(probU<plev) Idc<-"No  "
        else if(probL<plev&probU>=plev) Idc<-"?   "
        else if(probL>=plev) Idc<-"YifD"
	if(PFx>=PFx95U) Idc<-"Yes "
      }
      else if(Id==1) { # type-1 changepoints
        if(PFx<PFx95L) Idc<-"No  "
        else if(PFx>=PFx95L&PFx<PFx95U) Idc<-"?   "
        else if(PFx>=PFx95U) Idc<-"Yes "
      }

      cat(paste(sprintf("%1.0f",as.numeric(Id))," ",
              sprintf("%-4.4s",Idc),
              sprintf("%10.0f",IY0[Ic])," (",
              sprintf("%6.4f",probL),"-",
              sprintf("%6.4f",probU),")",
	      sprintf("%6.3f",plev),
	      sprintf("%10.4f",PFx)," (",
	      sprintf("%10.4f",PFx95L),"-",
	      sprintf("%10.4f",PFx95U),")\n",sep=""),
	      file=ofileIout,
	      append=TRUE)
      cat(paste("PMF : c=", sprintf("%4.0f",Ic), 
    	      "; (Time ", sprintf("%10.0f",IY0[Ic]), 
	      "); Type= ",sprintf("%4.0f",as.numeric(Id)),
              "; p=", sprintf("%10.4f",prob),
	      "(", sprintf("%10.4f",probL), 
	      "-", sprintf("%10.4f",probU), 
	      "); PFmax=", sprintf("%10.4f",PFx), 
	      "; CV95=",sprintf("%10.4f",PFx95), 
	      "(", sprintf("%10.4f",PFx95L), 
	      "-", sprintf("%10.4f",PFx95U),
	      "); Nseg=", sprintf("%4.0f",Nseg),"\n",sep=""), 
	      file=ofileSout, append=T)
    }
  }
  if(GUI)
    return(0)
  else{
    file.copy(from=ofileIout,to=ofileMout,overwrite=TRUE)
    cat("FindUD finished successfully...\n")
  }
}

StepSize<-function(InSeries,InCs,output,MissingValueCode,GUI=FALSE,p.lev=0.95,
                   Iadj=10000,Mq=10,Ny4a=0){
  ErrorMSG<-NA
  assign("ErrorMSG",ErrorMSG,envir=.GlobalEnv)
  if(!p.lev%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
      ErrorMSG<<-paste("FindU: input p.lev",p.lev,"error\n",
                 get("ErrorMSG",env=.GlobalEnv),"\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
  }
  plev<-p.lev
  if(Ny4a>0&Ny4a<=5) Ny4a<-5
  pkth<-match(p.lev,c(0.75,0.8,0.9,0.95,0.99,0.9999))
  flog<-paste(output,".log",sep="")
  itmp<-Read(InSeries,MissingValueCode)
  if(itmp<0){
    ErrorMSG<<-paste("StepSize: Error in read data from",InSeries,"\n",
               get("ErrorMSG",env=.GlobalEnv),"\n")
    if(!GUI) cat(ErrorMSG)
    return(-1)
  }
  N<-length(Y0); Nadj<-Ny4a*Nt
  readPFtable(N,pkth)
  itmp<-readLines(InCs)
  if(length(itmp)>=2){
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
      ErrorMSG<<-paste("StepSize: Ips read in from ",InCs,"error:")
      for(i in 1:Ns)
        ErrorMSG<<-paste(get("ErrorMSG",env=.GlobalEnv),Ips[i])
      ErrorMSG<<-paste(get("ErrorMSG",env=.GlobalEnv),"\n\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
    }
  }
  else{
    Ips<-N
    Ns<-0
  }
  Pk0<-Pk.PMFT(N)
  oout<-rmCycle(itable)
  Y1<-oout$Base
  EB<-oout$EB
  assign("EB",EB,envir=.GlobalEnv)
  if(length(EB)!=length(Icy)) {
    ErrorMSG<<-paste("Annual cycle length (from non-missing data) differ from original dataset",
                     "\n",get("ErrorMSG",env=.GlobalEnv),"\n")
    if(!GUI) cat(ErrorMSG)
    return(-1)
  }

  if(Ns>0) {
    Nsegs<-Ips-c(0,Ips[1:Ns])
    Iseg.longest<-sort(Nsegs,index=T,decreasing=T)$ix[1]
  } 
  else Iseg.longest<-0

  if(Iadj>(Ns+1)|Iseg.longest==0) Iseg.adj<-Ns+1 
  else if(Iadj==0)Iseg.adj<-Iseg.longest
  else Iseg.adj<-Iadj

  Ip0<-N
  otmp<-LSmultiRedCycle(Y1,Ti,Ip0,1)
  beta0<-otmp$trend
  betaL0<-otmp$betaL
  betaU0<-otmp$betaU
  meanhat0<-otmp$meanhat
  Ehat0<-mean(meanhat0)
  corD<-otmp$cor
  corDL<-otmp$corl
  corDU<-otmp$corh
  p.tr0<-otmp$p.tr

  otmp<-LSmultiRedCycle(Y1,Ti,Ips,Iseg.adj)
  Y1<-otmp$Y0
  cor<-otmp$cor
  corl<-otmp$corl
  corh<-otmp$corh
  pcor<-pt(abs(cor)*sqrt((N-2)/(1-cor^2)),N-2)
  Rf<-otmp$resi
  W<-otmp$W
  WL<-otmp$WL
  WU<-otmp$WU
  EB1<-otmp$EB

  itmp1<-cbind(EB1,Icy)
  itmp2<-cbind(1:N,Imd)
  colnames(itmp2)<-c("idx","Icy")
  itmp<-merge(itmp1,itmp2,by="Icy")
  EBfull<-itmp[order(itmp[,"idx"]),"EB1"]
  EEB<-mean(EB1,na.rm=T)

  if(Ns>0){
    Rb<-Y1-otmp$trend*Ti+EBfull
    QMout<-QMadjGaussian(Rb,Ips,Mq,Iseg.adj,Nadj)
    B<-QMout$PA
  }
  else B<-Y1-otmp$trend*Ti+EBfull

  adj<-otmp$Y0+EBfull
  B<-B+otmp$trend*Ti

  ofileAout<-paste(output,"_F.dat",sep="")
  ofilePdf<-paste(output,"_F.pdf",sep="")
  ofileSout<-paste(output,"_Fstat.txt",sep="")
  ofileIout<-paste(output,"_fCs.txt",sep="")
  ofileMout<-paste(output,"_mCs.txt",sep="")
  file.create(ofileSout)
  file.create(ofileAout)
  file.create(ofileIout)
  file.create(ofilePdf)
  cat(paste("The nominal level of confidence (1-alpha)=",plev,"\n"),file=ofileSout)
  cat(paste("Input data filename:", InSeries,"N=",N,"\n"),file=ofileSout,append=T)
  cat(file=ofileSout,paste(" Ignore changepoints -> trend0 =",
      round(beta0,6),"(",round(betaL0,6),",",round(betaU0,6),
      ") (p=",round(p.tr0,4),"); cor=",round(corD,4),"(",round(corDL,4),
      ",",round(corDU,4),")\n"),append=T)
  cat("Common trend TPR fit to the seasonal-cycle-adjusted Base series:\n",
      file=ofileSout,append=T)
  if(Ns>0) cat(paste("Nseg_shortest =",QMout$Nseg.mn,"; Mq = ",QMout$Mq,"; Ny4a = ",Ny4a,"\n"),
      file=ofileSout,append=T)
  if(Ns>0){
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
        cat(paste("Seg. ",i,": mean of QM-adjustments =",round(mean(QMout$W[I1:I2]),4),
            "\n",sep=""),file=ofileSout,append=T)
      }
    }
  }
  cat(paste("#steps= ",Ns,"; trend=",round(otmp$trend,6),"(",
            round(otmp$betaL,6),",",round(otmp$betaU,6),") (p=",
            round(otmp$p.tr,4),"); cor=",
            round(cor,4),"(",round(corl,4),",",round(corh,4),")",
	    round(pcor,4),"\n"),
	    file=ofileSout,append=T)
  if(Ns>0) for(i in 1:(Ns+1)){
    I1<-if(i==1) 1 else Ips[i-1]+1
    I2<-Ips[i]
    Delta<-otmp$mu[Iseg.adj]-otmp$mu[i]
    adj[I1:I2]<-adj[I1:I2]+Delta
    stepsize<-otmp$mu[i+1]-otmp$mu[i]
    cat(paste(Ips[i],IY0[Ips[i]],"stepsize=",round(stepsize,4),"\n"),
        file=ofileSout,append=T)
  }

  oR<-Y1-otmp$meanhat
  oR[2:N]<-oR[2:N]-oR[1:(N-1)]*cor
  Ehat<-mean(otmp$meanhat)
  meanhat0<-meanhat0-Ehat0+Ehat

  pdf(file=ofilePdf,onefile=T,paper='letter')
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  par(mfrow=c(2,1),cex.main=.8,cex.lab=.8,cex.axis=.8,cex=.8)

  uyrs<-unique(floor(ori.itable[,1]/10))*10
  labels<-NULL
  ats<-NULL
  for(i in 1:length(uyrs)){
    if(!is.na(match(uyrs[i],ori.itable[,1]))){
      labels<-c(labels,uyrs[i])
      ats<-c(ats,match(uyrs[i],ori.itable[,1]))
    }
  }
  par(mar=c(3,4,3,2)+.1)
  pdata<-rep(NA,dim(ori.itable)[1])
  pdata[ooflg]<-otmp$Y0
  plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(otmp$Y0,otmp$meanhat),max(otmp$Y0,otmp$meanhat)),
       xaxt="n",col="black",lwd=0.5,
       main="Base anomaly series and regression fit")
  axis(side=1,at=ats,labels=labels)
  pdata[ooflg]<-otmp$meanhat
  lines(1:dim(ori.itable)[1],pdata,col="red")

  pdata[ooflg]<-otmp$Y0+EBfull
  plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(otmp$Y0+EBfull,otmp$meanhat+EBfull),
              max(otmp$Y0+EBfull,otmp$meanhat+EBfull)),
       xaxt="n",col="black",lwd=0.5,
       main="Base series and regression fit")
  axis(side=1,at=ats,labels=labels)

  pdata[ooflg]<-otmp$meanhat+EEB
  lines(1:dim(ori.itable)[1],pdata,col="red")
  
  pdata[ooflg]<-adj
  plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(c(adj,B)),max(c(adj,B))),
       xaxt="n",col="black",lwd=0.5,
       main="Mean-adjusted base series")
  axis(side=1,at=ats,labels=labels)

  if(Ns>0) if(QMout$Mq>1){
    pdata[ooflg]<-B
    plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
         ylim=c(min(c(adj,B)),max(c(adj,B))),
         xaxt="n",col="black",lwd=0.5,
         main="QM-adjusted base series")
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

  odata<-matrix(NA,dim(ori.itable)[1],10)
# odata[ooflg,1]<-Ti
  odata[,1]<-c(1:dim(ori.itable)[1])
  odata[,2]<-ori.itable[,1]*10000+ori.itable[,2]*100+ori.itable[,3]
# odata[ooflg,3]<-round(otmp$Y0+EBfull,4)
  odata[,3]<-ori.itable[,4]
  odata[ooflg,4]<-round(otmp$meanhat+EEB,4)
  odata[ooflg,5]<-round(adj,4)
  odata[ooflg,6]<-round(otmp$Y0,4)
  odata[ooflg,7]<-round(otmp$meanhat,4)
  odata[ooflg,8]<-round(otmp$meanhat+EBfull,4)
# odata[ooflg,9]<-round(oR,4)
  if(Ns>0) if(QMout$Mq>1) odata[ooflg,9]<-round(B,4)
  odata[ooflg,10]<-round(meanhat0,4)

  Imd1<-ori.itable[,2]*100+ori.itable[,3]
  if(sum(is.na(ori.itable[,4])==F&Imd1==229)>0){
    if(Ns>0){
      tdata<-ori.itable[is.na(ori.itable[,4])==F,]
      IY1<-tdata[,1]*10000+tdata[,2]*100+tdata[,3]
      Ips.ymd<-IY0[Ips]
      Ips.1<-rep(NA,Ns+1)
      for(i in 1:Ns) Ips.1[i]<-c(1:length(IY1))[IY1==Ips.ymd[i]]
      Ips.1[Ns+1]<-length(IY1)
#     Ips.1<-c(1:length(IY1))[Ips.ymd==IY1]
      Imd2<-tdata[,2]*100+tdata[,3]
      Ids.leap<-c(1:length(Imd2))[Imd2==229]
      Nl<-length(Ids.leap)
      Rb<-Y1-otmp$trend*Ti+EBfull
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
        odata[is.na(odata[,3])==F&Imd1==229,9]<-round(B1.leap,4)
      }
    }
    else
      odata[Imd1==229,9]<-odata[Imd1==229,3]
    Ids.leapo<-c(1:dim(ori.itable)[1])[is.na(ori.itable[,4])==F&Imd1==229]
    for(jth in 1:length(Ids.leapo)){
      kth<-Ids.leapo[jth]
      if(Ns>0){
        k1th<-if(odata[kth-1,2]%in%IY0[Ips]) (kth+1) else (kth-1)
      }
      else k1th<-kth-1
      for(pth in c(4,7,8,10)) odata[kth,pth]<-odata[k1th,pth]
      for(pth in c(5,6)){delta1<-odata[k1th,3]-odata[k1th,pth]; odata[kth,pth]<-odata[kth,3]-delta1}
    }
  }
    
  write.table(file=ofileAout,odata,na=MissingValueCode,
	col.names=F,row.names=F)
  otmp<-LSmultiple(W,Ti,Ips)
  RW<-otmp$resi
  otmp<-LSmultiple(WL,Ti,Ips)
  RWL<-otmp$resi
  otmp<-LSmultiple(WU,Ti,Ips)
  RWU<-otmp$resi
  if(Ns==0) {
    cat(paste(Ns,"changepoints in Series", InSeries,"\n"),
        file=ofileIout)
#   ErrorMSG<<-"StepSize: PMFT finds the series to be homogeneous!\n\n"
#   if(!GUI) cat(ErrorMSG)
#   return(-2)
  }
  else{
    cat(paste(Ns,"changepoints in Series", InSeries,"\n"),
        file=ofileIout)

    for(i in 1:Ns){
      I1<- if(i==1) 1 else Ips[i-1]+1
      I3<-Ips[i+1]
      Ic<-Ips[i]
      Id<-Ids[i]
      Nseg<-I3-I1+1
    
      PFx95<-getPFx95(cor,Nseg)
      PFx95L<-getPFx95(corl,Nseg)
      PFx95U<-getPFx95(corh,Nseg)
      SSEf.Iseg<-sum(Rf[I1:I3]^2)
      Ips0<-Ips[-i]
      otmp<-LSmultiple(Y1,Ti,Ips0)
      SSE0.Iseg<-sum(otmp$resi[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      Pk0<-Pk.PMFT(Nseg)
      PFx<-Fx*Pk0[Ic-I1+1]

      otmp<-LSmultiple(W,Ti,Ips0)
      SSE0.Iseg<-sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg<-sum(RW[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(Fx<0){
        PFx<-0
        Fx<-0
        prob<-0
      }
      else prob<-pf(Fx,1,Nseg-3)

      otmp<-LSmultiple(WL,Ti,Ips0)
      SSE0.Iseg<-sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg<-sum(RWL[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      probL0<-if(Fx<0) 0 else pf(Fx,1,Nseg-3)

      otmp<-LSmultiple(WU,Ti,Ips0)
      SSE0.Iseg<-sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg<-sum(RWU[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      probU0<-if(Fx<0) 0 else pf(Fx,1,Nseg-3)

      probL<-min(probL0,probU0)
      probU<-max(probL0,probU0)

      if(Id==0) { # type-0 changepoints
        if(probU<plev) Idc<-"No  "
        else if(probL<plev&probU>=plev) Idc<-"?   "
        else if(probL>=plev) Idc<-"YifD"
	if(PFx>=PFx95U) Idc<-"Yes "
      }
      else if(Id==1) { # type-1 changepoints
        if(PFx<PFx95L) Idc<-"No  "
        else if(PFx>=PFx95L&PFx<PFx95U) Idc<-"?   "
        else if(PFx>=PFx95U) Idc<-"Yes "
      }

      cat(paste(sprintf("%1.0f",Id)," ",
                sprintf("%-4.4s",Idc),
                sprintf("%10.0f",IY0[Ic])," (",
	        sprintf("%10.4f",probL),"-",
	        sprintf("%10.4f",probU),")",
		sprintf("%6.3f",plev),
	        sprintf("%10.4f",PFx)," (",
	        sprintf("%10.4f",PFx95L),"-",
	        sprintf("%10.4f",PFx95U),")\n",sep=""),file=ofileIout,
	        append=TRUE)
      cat(paste("PMF : c=", sprintf("%4.0f",Ic), 
    	        "; (Time ", sprintf("%10.0f",IY0[Ic]),
	        "); Type= ",sprintf("%4.0f",Id),
                "; p=", sprintf("%10.4f",prob), 
	        "(", sprintf("%10.4f",probL), 
	        "-", sprintf("%10.4f",probU), 
	        "); PFmax=", sprintf("%10.4f",PFx), 
	        "; CV95=", sprintf("%10.4f",PFx95), 
	        "(", sprintf("%10.4f",PFx95L), 
	        "-", sprintf("%10.4f",PFx95U),
	        "); Nseg=", sprintf("%4.0f",Nseg),
	        "\n",sep=""), file=ofileSout, append=T)
    }
  }
  if(GUI)
    return(0)
  else{
    file.copy(from=ofileIout,to=ofileMout,overwrite=TRUE)
    cat("StepSize finished successfully...\n")
  }
}

QMadj.GaussianDLY<-function(InSeries,InCs,output,MissingValueCode,GUI=FALSE,
    Iadj=10000,Mq=10,Ny4a=0){
  ErrorMSG<-NA
  assign("ErrorMSG",ErrorMSG,envir=.GlobalEnv)
  itmp<-ReadDLY.g(InSeries,MissingValueCode)
  if(itmp<0){
    ErrorMSG<<-paste("QMadj.GaussianDLY: Error in read data from",InSeries,"\n",
               get("ErrorMSG",env=.GlobalEnv),"\n")
    if(!GUI) cat(ErrorMSG)
    return(-1)
  }
  N<-length(Y0); Nadj<-Ny4a*Nt
  itmp<-readLines(InCs)
  if(length(itmp)>=2){
    Ns<-length(itmp)-1
    Ips<-c(rep(0,Ns),N)
    for(i in 1:Ns){ # using YYYYMMDD as index, searching for the largest
                    # date less or equal to given YYYYMMDD
      ymdtmp<-as.numeric(substr(itmp[i+1],7,16))
      if(ymdtmp==as.integer(ymdtmp/100)*100) ymdtmp<-ymdtmp+15
      # set date as 15 if input break point is 00 for date
      it<-match(ymdtmp,IY0)
      if(!is.na(it)) Ips[i]<-it
      else Ips[i]<-max(c(1:N)[IY0<=ymdtmp])
    }
    if(sum(is.na(Ips))>0|!identical(Ips,sort(Ips))){
      ErrorMSG<<-paste("QMadj.GaussianDLY: Ips read in from ",InCs,"error:")
      for(i in 1:Ns)
        ErrorMSG<<-paste(get("ErrorMSG",env=.GlobalEnv),Ips[i])
      ErrorMSG<<-paste(get("ErrorMSG",env=.GlobalEnv),"\n\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
    }
  }
  else Ns<-0

  ofileSout<-paste(output,"_QMadjDLYstat.txt",sep="")
  file.create(ofileSout)
  cat(paste("Input data filename:", InSeries,"; N=",N,"\n"),file=ofileSout)

  if(Ns>0) {
    Nsegs<-Ips-c(0,Ips[1:Ns])
    Iseg.longest<-sort(Nsegs,index=T,decreasing=T)$ix[1]
  }
  else Iseg.longest<-0

  if(Iadj>(Ns+1)|Iseg.longest==0) Iseg.adj<-Ns+1
  else if(Iadj==0)Iseg.adj<-Iseg.longest
  else Iseg.adj<-Iadj

  oout<-rmCycle(itable)
  Y1<-oout$Base
  EB<-oout$EB
  assign("EB",EB,envir=.GlobalEnv)
  if(length(EB)!=length(Icy)) {
    ErrorMSG<<-paste("Annual cycle length (from non-missing data) differ from original dataset",
                     "\n",get("ErrorMSG",env=.GlobalEnv),"\n")
    if(!GUI) cat(ErrorMSG)
    return(-1)
  }

  otmp<-LSmultiRedCycle(Y1,Ti,Ips,Iseg.adj)
  Y1<-otmp$Y0
  cor<-otmp$cor
  corl<-otmp$corl
  corh<-otmp$corh
  pcor<-pt(abs(cor)*sqrt((N-2)/(1-cor^2)),N-2)
  Rf<-otmp$resi
  W<-otmp$W
  L<-otmp$WL
  WU<-otmp$WU
  EB1<-otmp$EB

  itmp1<-cbind(EB1,Icy)
  itmp2<-cbind(1:N,Imd)
  colnames(itmp2)<-c("idx","Icy")
  itmp<-merge(itmp1,itmp2,by="Icy")
  EBfull<-itmp[order(itmp[,"idx"]),"EB1"]
  EEB<-mean(EB1,na.rm=T)

  if(Ns>0){
    Rb<-Y1-otmp$trend*Ti+EBfull
#   QMout<-QMadjDLY(Rb,Ips,Mq,Iseg.adj)
    QMout<-QMadjGaussian(Rb,Ips,Mq,Iseg.adj,Nadj)
    B<-QMout$PA
  }
  else B<-Y1-otmp$trend*Ti+EBfull

  adj<-otmp$Y0+EBfull
  B<-B+otmp$trend*Ti

  if(Ns>0){
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
        cat(paste("Seg. ",i,": mean of QM-adjustments =",round(mean(QMout$W[I1:I2]),4),
            "\n",sep=""),file=ofileSout,append=T)
      }
    }
  }

  cat(paste("#steps= ",Ns,"; trend=",round(otmp$trend,6),"(",
            round(otmp$betaL,6),",",round(otmp$betaU,6),") (p=",
	    round(otmp$p.tr,4),"); cor=",
	    round(cor,4),"(",round(corl,4),",",round(corh,4),")",
	    round(pcor,4),"\n"),
	    file=ofileSout,append=T)

  if(Ns>0) for(i in 1:(Ns+1)){
    I1<-if(i==1) 1 else Ips[i-1]+1
    I2<-Ips[i]
    Delta<-otmp$mu[Iseg.adj]-otmp$mu[i]
    adj[I1:I2]<-adj[I1:I2]+Delta
    stepsize<-otmp$mu[i+1]-otmp$mu[i]
    cat(paste(Ips[i],IY0[Ips[i]],"stepsize=",round(stepsize,4),"\n"),
        file=ofileSout,append=T)
  }

# oR<-Y1-otmp$meanhat
# oR[2:N]<-oR[2:N]-oR[1:(N-1)]*cor
# Ehat<-mean(otmp$meanhat)

  if(Ns>0){
    odata<-matrix(NA,dim(ori.itable)[1],8)
    odata[,1]<-ori.itable[,1]
    odata[,2]<-ori.itable[,2]
    odata[,3]<-ori.itable[,3]
    odata[olflg,4]<-round(otmp$Y0+EBfull,4)
    if(QMout$Mq>1) odata[olflg,5]<-round(B,4)
    odata[olflg,6]<-round(adj,4)
    if(QMout$Mq>1) odata[olflg,7]<-round(B-otmp$Y0-EBfull,4)
    odata[olflg,8]<-round(adj-otmp$Y0-EBfull,4)
  }
  else odata<-cbind(ori.itable[,c(1,2,3,4,4,4)],0,0)

  ofileAout<-paste(output,"_QMadjDLY.dat",sep="")
# ofileAout<-output # suggested by Lucie, user can choose whatever filename
  write.table(odata,file=ofileAout,na=MissingValueCode,col.names=F,row.names=F)

  ofilePdf<-paste(output,"_QMadjDLY.pdf",sep="")
  pdf(file=ofilePdf,onefile=T,paper='letter')
  op <- par(no.readonly = TRUE) # the whole list of settable par's
  par(mfrow=c(2,1),cex.main=.8,cex.lab=.8,cex.axis=.8,cex=.8)

  uyrs<-unique(floor(ori.itable[,1]/10))*10
  labels<-NULL
  ats<-NULL
  for(i in 1:length(uyrs)){
    if(!is.na(match(uyrs[i],ori.itable[,1]))){
      labels<-c(labels,uyrs[i])
      ats<-c(ats,match(uyrs[i],ori.itable[,1]))
    }
  }
  par(mar=c(3,4,3,2)+.1)
  pdata<-rep(NA,dim(ori.itable)[1])
  pdata[olflg]<-otmp$Y0
  plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(otmp$Y0,otmp$meanhat),max(otmp$Y0,otmp$meanhat)),
       xaxt="n",col="black",lwd=.3,
       main="Base anomaly series and regression fit")
  axis(side=1,at=ats,labels=labels)
  pdata[olflg]<-otmp$meanhat
  lines(1:dim(ori.itable)[1],pdata,col="red")

  pdata[olflg]<-otmp$Y0+EBfull
  plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(otmp$Y0+EBfull,otmp$meanhat+EBfull),
              max(otmp$Y0+EBfull,otmp$meanhat+EBfull)),
       xaxt="n",col="black",lwd=.3,
       main="Base series and regression fit")
  axis(side=1,at=ats,labels=labels)

  pdata[olflg]<-otmp$meanhat+EEB
  lines(1:dim(ori.itable)[1],pdata,col="red")
  
  pdata[olflg]<-adj
  plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(c(adj,B)),max(c(adj,B))),
       xaxt="n",col="black",lwd=.3,
       main="Mean-adjusted base series")
  axis(side=1,at=ats,labels=labels)

  if(Ns>0) if(QMout$Mq>1){
    pdata[olflg]<-B
    plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
         ylim=c(min(c(adj,B)),max(c(adj,B))),
         xaxt="n",col="black",lwd=.3,
         main="QM-adjusted base series")
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
}

ReadDLY.g<-function(idata,MissingValue){
  if(!file.exists(idata)) {
    ErrorMSG<<-paste("Input datafile",idata,"does not exist!\n")
    return(-1)
  }
  if(is.csv(idata)){
    itmp<-try(read.table(idata,sep=",",header=F,na.strings=MissingValue,
            colClasses=rep("real",4)),silent=T)
    if(inherits(itmp,"try-error")){
      ErrorMSG<<-geterrmessage()
      return(-1)
    }
    else itable<-itmp
  }
  else{
    itmp<-try(read.table(idata,sep="",header=F,na.strings=MissingValue,
            colClasses=rep("real",4)),silent=T)
    if(inherits(itmp,"try-error")){
      ErrorMSG<<-geterrmessage()
      return(-1)
    }
    else itable<-itmp
  }
  if(ncol(itable)!=4){
    ErrorMSG<<-paste(idata,"has",ncol(itable),"columns. The number of columns should be 4\n")
    return(-1)
  }
  colnames(itable)<-c("id1","id2","id3","data")

  iyrbegin<-itable[1,1]
  imdbegin<-itable[1,2]*100+itable[1,3]
  iyrend<-itable[dim(itable)[1],1]
  imdend<-itable[dim(itable)[1],2]*100+itable[dim(itable)[1],3]
# keep input base data as ori.itable
  ori.itable<-itable
# check input data (both base and ref), no jump with begin and end
  Icy<-sort(unique(itable[,2]*100+itable[,3]))
  Ind2<-iyrbegin*10000+Icy[Icy>=imdbegin] # first year
# if(iyrend>(iyrbegin+1)) for(i in (iyrbegin+1):(iyrend-1))
#   Ind2<-c(Ind2,i*10000+Icy)
# Ind2<-c(Ind2,iyrend*10000+Icy[Icy<=imdend])
  Nt<-length(Icy)
  Nall<-dim(itable)[1]
  ind<-itable[,1]*10000+itable[,2]*100+itable[,3]
# for(i in 1:length(Ind2)) 
#   if(Ind2[i]!=ind[i]) {
#     ErrorMSG<<-paste("Input data series not continuous at:",Ind2[i],ind[i],"\n")
#     return(-1)
#   }
  IY0<-ind[is.na(itable[,4])==F]
  IY0flg<-rep(0,length(IY0))
  Y0<-itable[is.na(itable[,4])==F,4]
  olflg<-is.na(itable[,4])==F
  Iyr<-floor(IY0/10000)
  Imd<-IY0-Iyr*10000
  Ti<-IY0
  for(i in 1:length(IY0)){
    ith<-match(Imd[i],Icy)
    Ti[i]<-(Iyr[i]-iyrbegin)*Nt+ith
  }
  for(i in 1:(length(IY0)-1)){
    if(Ti[i+1]-Ti[i]==1) IY0flg[i]<-1
  }
  if(sum(IY0flg)<10){  # too few available data for autocorlh 
    ErrorMSG<<-paste("Too many missing values in ", idata, "to estimate autocorrelation\n")
    return(-1)
  }
  itable<-itable[is.na(itable[,4])==F,]
  assign("ori.itable",ori.itable,envir=.GlobalEnv)
  assign("itable",itable,envir=.GlobalEnv)
  assign("Ti",Ti,envir=.GlobalEnv) # Time index for LS fitting
  assign("Y0",Y0,envir=.GlobalEnv) # Data series for Base
  assign("IY0",IY0,envir=.GlobalEnv) # Cycle index for Base
  assign("Imd",Imd,envir=.GlobalEnv) # Cycle index for Base
  assign("IY0flg",IY0flg,envir=.GlobalEnv) # continuous flag for Base
  assign("Icy",Icy,envir=.GlobalEnv) # Cycle index
  assign("Nt",Nt,envir=.GlobalEnv) # Cycle length
  assign("olflg",olflg,envir=.GlobalEnv)
  return(0)
}

QMadjGaussian<-function(P,Ips,Mq,Iseg.adj,Nadj){
  Ns<-length(Ips)-1
  N<-length(P)
  Nseg.mn<-N
  for(i in 1:(Ns+1)){
    I1<-if(i==1) 1 else Ips[i-1]+1
    I2<-Ips[i]
    Nseg<-I2-I1+1
    if(Nseg<Nseg.mn) Nseg.mn<-Nseg
  }
  if(Nt<=12) Ln<-5 else Ln<-20
  if(Mq<=0) Mq<-min(floor(Nseg.mn/Ln),100) else Mq<-min(floor(Nseg.mn/Ln),Mq)
  if(Mq>100) Mq<-100
  if(Mq<=0) Mq<-1
  Fd<-.5/Mq
  Fcat<-matrix(NA,(Ns+1),(Mq+2))
  F<-matrix(NA,(Ns+1),N)
  EPa<-matrix(0,(Ns+1),(Mq+2))
  EPb<-matrix(0,(Ns+1),(Mq+2))
  for(i in 1:(Ns+1)) Fcat[i,]<-seq(0,by=1/Mq,length=(Mq+2))
  for(i in 1:Ns){
    I1<-if(i==1) 0 else Ips[i-1]
    I2<-Ips[i]
    if(Nadj>0) I1<-max(c(I1,I2-Nadj))
    Nseg<-I2-I1
    Y<-P[(I1+1):I2]
    iindx<-sort(Y,index=T)$ix
    irank<-sort(iindx,index=T)$ix
#   F[i,1:Nseg]<-irank/Nseg
    EPb[i,]<-0
    for(k in 2:(Mq+1)){
      Mp1<-floor(Nseg*Fcat[i,(k-1)])
      Mp2<-floor(Nseg*Fcat[i,k])
      EPb[i,k]<-mean(Y[iindx[(Mp1+1):Mp2]])
    }
    EPb[i,(Mq+2)]<-EPb[i,(Mq+1)]
    I3<-Ips[i]
    I4<-Ips[i+1]
    if(Nadj>0) I4<-min(c(I4,I3+Nadj))
    if(Nadj>0|i==Ns){
      Nseg<-I4-I3
      Y<-P[(I3+1):I4]
      iindx<-sort(Y,index=T)$ix
      irank<-sort(iindx,index=T)$ix
      EPa[i,]<-0
      for(k in 2:(Mq+1)){
        Mp1<-floor(Nseg*Fcat[i,(k-1)])
	Mp2<-floor(Nseg*Fcat[i,k])
	EPa[i,k]<-mean(Y[iindx[(Mp1+1):Mp2]])
      }
      EPa[i,(Mq+2)]<-EPa[i,(Mq+1)]
    }
    EPa[i,(Mq+2)]<-EPa[i,(Mq+1)]
  }
  EPb[,1]<-EPb[,2]
  EPa[,1]<-EPa[,2]

  if(Nadj==0) if(Ns>1) EPa[1:(Ns-1),]<-EPb[2:Ns,]
  Adj<-matrix(0,(Ns+1),(Mq+2))
  for(k in 1:(Mq+2)){
    if(Iseg.adj>1) for(i in 1:(Iseg.adj-1)) Adj[i,k]<-sum(EPa[i:(Iseg.adj-1),k])-sum(EPb[i:(Iseg.adj-1),k])
    if(Iseg.adj<=Ns) for(i in (Iseg.adj+1):(Ns+1)) Adj[i,k]<-sum(EPb[Iseg.adj:(i-1),k])-sum(EPa[Iseg.adj:(i-1),k])
  }

  osmean<-c(1:(Mq+2)) # for plot purpose
  for(Iseg in c(1:(Ns+1)))
    osmean<-cbind(osmean,Fcat[Iseg,]-Fd,Adj[Iseg,])
#   osmean<-cbind(osmean,Fcat[Iseg,]-Fd,EP[Iseg.adj,]-EP[Iseg,])
  # output osmean is a 2*(Ns+1)+1 by Mq+2 matrix

  PA<-P
  W<-rep(NA,N)
  osp<-NULL
  for(i in 1:(Ns+1)){
    I1<-if(i==1) 1 else Ips[i-1]+1
    I2<-Ips[i]
    Nseg<-I2-I1+1
    Y<-P[I1:I2]
    iindx<-sort(Y,index=T)$ix
    irank<-sort(iindx,index=T)$ix
    F[i,1:Nseg]<-irank/Nseg
    if(i==Iseg.adj) PA[I1:I2]<-P[I1:I2] else{
      dx<-Fcat[i,]-Fd
      fdx<-Adj[i,]
#     fdx<-EP[Iseg.adj,]-EP[i,]
      if(Mq==1) fdx[1]<-fdx[2]
      fdx2<-splineN(dx,fdx,2E30,2E30)
      for(j in I1:I2) W[j]<-splintN(dx,fdx,fdx2,F[i,(j-I1+1)])
      PA[I1:I2]<-P[I1:I2]+W[I1:I2]
    }

    Rs<-F[i,1:Nseg]
    ors<-sort(Rs,index=T)$ix
    osp<-rbind(osp,cbind(I1:I2,Rs[ors],W[I1:I2][ors]))
  }
  oout<-list()
  oout$PA<-PA
  oout$W<-W
  oout$Mq<-Mq
  oout$Nseg.mn<-Nseg.mn
  oout$osmean<-osmean
  oout$osp<-osp
  return(oout)
}

QMadjDLY<-function(P,Ips,Mq,Iseg.adj){
  Ns<-length(Ips)-1
  N<-length(P)
  Nseg.mn<-N
  for(i in 1:(Ns+1)){
    I1<-if(i==1) 1 else Ips[i-1]+1
    I2<-Ips[i]
    Nseg<-I2-I1+1
    if(Nseg<Nseg.mn) Nseg.mn<-Nseg
  }
  if(Mq<=0) Mq<-min(floor(Nseg.mn/20),100) else Mq<-min(floor(Nseg.mn/20),Mq)
  if(Mq>100) Mq<-100
  if(Mq<=0) Mq<-1
  Fd<-.5/Mq
  Fcat<-matrix(NA,(Ns+1),(Mq+2))
  F<-matrix(NA,(Ns+1),N)
  EP<-matrix(NA,(Ns+1),(Mq+2))
  for(i in 1:(Ns+1)) Fcat[i,]<-seq(0,by=1/Mq,length=(Mq+2))
  for(i in 1:(Ns+1)){
    I1<-if(i==1) 1 else Ips[i-1]+1
    I2<-Ips[i]
    Nseg<-I2-I1+1
    Y<-P[I1:I2]
    iindx<-sort(Y,index=T)$ix
    irank<-sort(iindx,index=T)$ix
    F[i,1:Nseg]<-irank/Nseg
    EP[i,]<-0
    for(k in 2:(Mq+1)){
      Mp1<-floor(Nseg*Fcat[i,(k-1)])
      Mp2<-floor(Nseg*Fcat[i,k])
      EP[i,k]<-mean(Y[iindx[(Mp1+1):Mp2]])
    }
    EP[i,(Mq+2)]<-EP[i,(Mq+1)]
  }
  EP[,1]<-EP[,2]

  osmean<-c(1:(Mq+2)) # for plot purpose
  for(Iseg in c(1:(Ns+1)))
    osmean<-cbind(osmean,Fcat[Iseg,]-Fd,EP[Iseg.adj,]-EP[Iseg,])
  # output osmean is a 2*(Ns+1)+1 by Mq+2 matrix

  PA<-P
  W<-rep(NA,N)
  osp<-NULL
  for(i in 1:(Ns+1)){
    I1<-if(i==1) 1 else Ips[i-1]+1
    I2<-Ips[i]
    Nseg<-I2-I1+1
    if(i==Iseg.adj) PA[I1:I2]<-P[I1:I2] else{
      dx<-Fcat[i,]-Fd
      fdx<-EP[Iseg.adj,]-EP[i,]
      if(Mq==1) fdx[1]<-fdx[2]
      fdx2<-splineN(dx,fdx,2E30,2E30)
      for(j in I1:I2) W[j]<-splintN(dx,fdx,fdx2,F[i,(j-I1+1)])
      PA[I1:I2]<-P[I1:I2]+W[I1:I2]
    }

    Rs<-F[i,1:Nseg]
    ors<-sort(Rs,index=T)$ix
    osp<-rbind(osp,cbind(I1:I2,Rs[ors],W[I1:I2][ors]))
  }
  oout<-list()
  oout$PA<-PA
  oout$W<-W
  oout$Mq<-Mq
  oout$Nseg.mn<-Nseg.mn
  oout$osmean<-osmean
  oout$osp<-osp
  return(oout)
}

QMadjGaussian.wRef<-function(B,Dif,Ips,Mq,Iseg.adj,Nadj,Nt,Ny4a){
  Ns<-length(Ips)-1
  N<-length(B)
  Nmin.Mq<- if(Nt<=12) 5 else 20
  
  Nseg.mn<-N

  for(i in 1:(Ns+1)){
    I1<-if(i==1) 1 else Ips[i-1]+1
    I2<-Ips[i]
    Nseg<-0
    DD<-Dif[I1:I2]
    Nseg<-sum(ifelse(is.na(DD), 0, 1))
    if(Nseg<Nseg.mn) Nseg.mn<-Nseg
  }

  if(Mq<0) { # not include equal
    Mq<-min(floor(Nseg.mn/Nmin.Mq),100) 
#   if(Nt<=12) Mq<-min(floor(Nseg.mn/5),20) 
    } else {
    Mq1<-min(floor(Nseg.mn/Nmin.Mq),Mq)
#   if(Nt<=12) Mq1<-min(floor(Nseg.mn/5),Mq) 
    Mq<-Mq1
  }
  if(Mq>100) Mq<-100
  if(Mq==0) Mq<-1

  tflg<-T
  while(tflg){
    tflg<-FALSE
    fMq<-floor(Mq)
    Fd<-.5/fMq

    EBa<-matrix(0,(Ns+1),(Mq+2))
    EBb<-matrix(0,(Ns+1),(Mq+2))

    Fcat<-matrix(0,(Ns+1),(Mq+2))
    for(i in 1:(Ns+1)) Fcat[i,]<-seq(0,by=1/Mq,length=(Mq+2))
  
    Nadj<-Ny4a*Nt
    Yd<-NULL
    F<-NULL
    EYd<-rep(0,Ns)

    for(i in 1:Ns){
      I1<-if(i==1) 0 else Ips[i-1]
      I2<-Ips[i]
      if(Ny4a>0) I1<-max(c(I1,I2-Nadj))
      Nseg<-I2-I1
      Y<-B[(I1+1):I2]
      Yd<-Dif[(I1+1):I2]
      EYd1<-mean(Yd,na.rm=TRUE)   
    
      iindx<-rev(sort(Y,index=T,decreas=T)$ix)
      irank<-sort(iindx,index=T)$ix
  
      for(k in 2:(Mq+1)){
        Mp1<-floor(Nseg)*Fcat[i,(k-1)] #? floor should ()the whole expression?
        Mp2<-floor(Nseg)*Fcat[i,k]
        if(abs(Mp1-round(Mp1))<1E-9) Mp1<-round(Mp1)
        if(abs(Mp2-round(Mp2))<1E-9) Mp2<-round(Mp2)
        Mp1<-floor(Mp1); Mp2<-floor(Mp2)
        if(sum(!is.na(Yd[iindx[(Mp1+1):Mp2]]))<Nmin.Mq) {
          if(Mq>1){
            tflg<-TRUE
            Mq<-Mq-1
          }
          else{
            ErrorMSG<<-paste("!!!! Warning - The Ref series has too many missing values to do QM-adjustments.\n",
            "Only mean-adjustments were done. You may want to use a better Ref series if possible!\n",
            get("ErrorMSG",env=.GlobalEnv),"\n")
          }
        }
        else EBb[i,k]<-mean(Yd[iindx[(Mp1+1):Mp2]],na.rm=TRUE)
      }
      EBb[i,1]<-EBb[i,2]
      EBb[i,(Mq+2)]<-EBb[i,(Mq+1)]
    
      I3<-Ips[i]
      I4<-Ips[i+1]
      if(Ny4a>0) I4<-min(c(I4,I3+Nadj))
      if(Iseg.adj>0|i>=Ns){
        Nseg<-I4-I3
        Y<-B[(I3+1):I4]
        Yd<-Dif[(I3+1):I4]
        EYd2<-mean(Yd,na.rm=TRUE)   
      }
      
      EYd[i]<-EYd2-EYd1
      
      iindx<-rev(sort(Y,index=T,decreas=T)$ix)
      irank<-sort(iindx,index=T)$ix
    
      for(k in 2:(Mq+1)){
        Mp1<-floor(Nseg)*Fcat[i,(k-1)]
        Mp2<-floor(Nseg)*Fcat[i,k]
        if(abs(Mp1-round(Mp1))<1E-9) Mp1<-round(Mp1)
        if(abs(Mp2-round(Mp2))<1E-9) Mp2<-round(Mp2)
        Mp1<-floor(Mp1); Mp2<-floor(Mp2)
        if(sum(!is.na(Yd[iindx[(Mp1+1):Mp2]]))<Nmin.Mq) {
          if(Mq>1){
            tflg<-TRUE
            Mq<-Mq-1
          }
          else{
            ErrorMSG<<-paste("!!!! Warning - The Ref series has too many missing values to do QM-adjustments.\n",
            "Only mean-adjustments were done. You may want to use a better Ref series if possible!\n",
            get("ErrorMSG",env=.GlobalEnv),"\n")
          }
        }
        EBa[i,k]<-mean(Yd[iindx[(Mp1+1):Mp2]],na.rm=TRUE)
      }
      EBa[i,1]<-EBa[i,2]
      EBa[i,(Mq+2)]<-EBa[i,(Mq+1)]
    }
  
    if(tflg==FALSE){
      AdjM<-rep(0,Ns+1)

      if(Iseg.adj>1) AdjM[1:(Iseg.adj-1)]<-rev(cumsum(EYd[(Iseg.adj-1):1]))
      if(Iseg.adj<=Ns) AdjM[(Iseg.adj+1):(Ns+1)]<- -cumsum(EYd[Iseg.adj:Ns])
   
      if(Ny4a==0 & Ns>1) EBa[1:(Ns-1),]<-EBb[2:Ns,]
 
      Adj<-matrix(0,Ns+1,(Mq+2))  

      for(k in 1:(Mq+2)){
        if(Iseg.adj>1) for(i in 1:(Iseg.adj-1)) Adj[i,k]<-sum(EBa[i:(Iseg.adj-1),k])-sum(EBb[i:(Iseg.adj-1),k])
        if(Iseg.adj<=Ns) for(i in (Iseg.adj+1):(Ns+1)) Adj[i,k]<-sum(EBb[Iseg.adj:(i-1),k])-sum(EBa[Iseg.adj:(i-1),k])
      }

      osmean<-c(1:(Mq+2)) 
      for(Iseg in c(1:(Ns+1))) osmean<-cbind(osmean,Fcat[Iseg,]-Fd,Adj[Iseg,])
  
      BA<-B
      osp<-NULL

      for(i in 1:(Ns+1)){
        I1<-if(i==1) 1 else Ips[i-1]+1
        I2<-Ips[i]
        Nseg<-I2-I1+1
    
        DB<-rep(NA,Nseg)
    
        if(i==Iseg.adj) BA[I1:I2]<-B[I1:I2] else{
          dx<-Fcat[i,]-Fd
          fdx<-Adj[i,]
          if(Mq==1) fdx[1]<-fdx[2]
          fdx2<-splineN(dx,fdx,2E30,2E30)
      
          iindx<-rev(sort(B[I1:I2],index=T,decreas=T)$ix)
          irank<-sort(iindx,index=T)$ix
      
          F<-rep(0,Nseg)
          F<-irank/Nseg
      
          dxi<-F
      
          for(k in 1:Nseg) DB[k]<-splintN(dx,fdx,fdx2,dxi[k])
     
          DBm0<-sum(DB)/Nseg
          remain<-AdjM[i]-DBm0
      
          DBpo<-DB[DB>0]
          DBne<-DB[DB<0]
       
          Npo<-length(DBpo)
          Nne<-length(DBne)
      
          sDBpo<-sum(DBpo)
          sDBne<-sum(DBne)
       
          if(remain<0) DB[DB<0]<-DB[DB<0]+remain*Nseg*DB[DB<0]/sDBne
          else DB[DB>0]<-DB[DB>0]+remain*Nseg*DB[DB>0]/sDBpo
           
          BA[I1:I2]<-B[I1:I2]+DB
          DBm<-sum(DB)/Nseg
      
          RO1<-AdjM[i]-DBm0
          RO2<-AdjM[i]-DBm
      
          if(i!=Iseg.adj){
            Rs<-F
            ors<-sort(Rs,index=T)$ix
            osp<-rbind(osp,cbind(I1:I2,Rs[ors],DB[ors]))
          }
        }
      }
    } # end if(!tflg)
  } # end while
  oout<-list()
  oout$PA<-BA
  oout$DB<-DB
  oout$Mq<-Mq
  oout$Nseg.mn<-Nseg.mn
  oout$osmean<-osmean
  oout$osp<-osp
  oout$AdjM<-AdjM
  return(oout)
}

splineN<-function(x,y,yp1,yp2){
  n<-length(x)
  if(length(y)!=n) stop("input vector length differ")
  y2<-rep(NA,0)
  u<-rep(NA,0)
  if(yp1>1E30){
    y2[1]<-0
    u[1]<-0 }
  else{
    y2[1]<-(-0.5)
    u[1]<-(3/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1)
  }
  if(n>2) for(i in c(2:(n-1))){
    sig<-(x[i]-x[i-1])/(x[i+1]-x[i-1])
    p<-sig*y2[i-1]+2
    y2[i]<-(sig-1)/p
    u[i]<-(6*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))/
          (x[i+1]-x[i-1])-sig*u[i-1])/p
  }
  if(yp2>1E30){
    qn<-0
    un<-0 }
  else{
    qn<-0.5
    un<-(3/(x[n]-x[n-1]))*(yp2-(y[n]-y[n-1])/(x[n]-x[n-1]))
  }
  y2[n]<-(un-qn*u[n-1])/(qn*y2[n-1]+1)
  for(i in c((n-1):1)) y2[i]<-y2[i]*y2[i+1]+u[i]
  return(y2)
}

splintN<-function(xa,ya,y2a,x){
  n<-length(xa)
  if(length(ya)!=n|length(y2a)!=n) stop("input vector length differ")
  klo<-1
  khi<-n
  while((khi-klo)>1){
    k<-ceiling((khi+klo)/2)
    if(xa[k]>x) khi<-k else klo<-k
  }
  h<-xa[khi]-xa[klo]
  if(h==0) stop("bad xa input in splintN")
  a<-(xa[khi]-x)/h
  b<-(x-xa[klo])/h
  y<-a*ya[klo]+b*ya[khi]+((a**3-a)*y2a[klo]+(b**3-b)*y2a[khi])*(h**2)/6
  return(y)
}



Pk.PMFT<-function(N){
# if(floor(N)!=N) stop("input data error in Pk")
  Nlt40<- if(N<40) TRUE else FALSE
  Nle500<- if(N<=500) TRUE else FALSE
  
  K<-seq(1,(N-1))
  Kmin<- if(floor((N-1)/2)==(N-1)/2) c(1:floor((N-1)/2),floor((N-1)/2):1)
         else c(1:(floor((N-1)/2)+1),floor((N-1)/2):1)
  W<- floor(if(Nle500) 11*N/50 else 21*N/100)
  A<-abs(1-2*K/N)
  B<-log(N)
  C<-log(B)
  D<-log(log(N+150))
  Q<-c(abs(1-50*K[1:floor(N/2)]/(N*11)),
       abs(1-(50*K[(floor(N/2)+1):(N-1)]-28*N)/(N*11)))
  tmp1<-11/B^3
  tmp2<-abs(2*C^4/(900-9*Kmin))
  flg<-tmp1<tmp2
  S<-tmp1*flg+tmp2*(!flg)
  S[Kmin==100]<-tmp1
  tmp1<-(1-Q^(B*(3-C)/6))
  tmp2<-(1+A^(B*(3-C)/6))
  F<-c(tmp1[1:W],tmp2[(W+1):(N-W-1)],tmp1[(N-W):(N-1)])
  tmp1<-(11.7-S)*B^.01-11.8
  tmp2<-((C+353)*N^.01-340)/200
  v<-c(tmp1[1:W],rep(tmp2,(N-W*2-1)),tmp1[(N-W):(N-1)])

  tmp1<-((64*N^(1/40)+35-C)*F^v)/100
  tmp2<-min(c(1,(98000+6*N)/100000))*F^v
  P0<-c(tmp1[1:W],tmp2[(W+1):(N-W-1)],tmp1[(N-W):(N-1)])

  P<-P0
  if(N>=40){
    L<-floor(1+(317*N^.75-2867)/1000)
    if(N<=100)
      delta<-rep(D^(1/3)*(P0[L+1]-P0[L])+C^3/(N*5),N-1)
    else{
      delta<-rep(NA,N-1)
      delta[1:L]<-(P0[L]-P0[1])*A[1:L]^(C^3)/(L+B-2*C-1)
      delta[(N-L):(N-1)]<-(P0[N-L]-P0[N-1])*A[(N-L):(N-1)]^(C^3)/(L+B-2*C-1)
    }
    P[1:L]<-P0[L]-(L-(1:L))*delta[1:L]
    P[(N-L):(N-1)]<-P0[N-L]-((N-L):(N-1)-N+L)*delta[(N-L):(N-1)]
  }
  return(P)
}

LSmatrix<-function(Y,T,Ic){
  Nx<-length(Y)
  D<-rep(1,Nx)
  X<-t(t(Y))
  D<-cbind(D,T-mean(T))
  if(!is.na(Ic)) D<-cbind(D,c(rep(0,Ic),rep(1,Nx-Ic)))
  sig<-solve(t(D)%*%D)%*%t(D)%*%X
  fitted<-D%*%sig
  resi<-X-fitted
  SSE<-sum(resi^2)
  oout<-list()
  oout$sig<-as.vector(sig)
  oout$fitted<-as.vector(fitted)
  oout$resi<-as.vector(resi)
  oout$SSE<-SSE
  return(oout)
}

PMFT<-function(Y,T,Pk0){
  N<-length(Y)
  PFx<-(-99999.)
  Fx<-(-99999.)
  KPx<-0
  oout1<-LSmatrix(Y,T,NA)
  for(i in Nmin:(N-Nmin)){
    oout2<-LSmatrix(Y,T,i)
    Fc<-(oout1$SSE-oout2$SSE)*(N-3)/oout2$SSE
    PFc<-Fc*Pk0[i]
    if(PFc>PFx){
      PFx<-PFc
      KPx<-i
      Fx<-Fc
    }
  }
  oout<-list()
  oout$PFx<-PFx
  oout$KPx<-KPx
  oout$Fx<-Fx
  return(oout)
}

PMFxKc<-function(Y,T,I0,I2,Ic){
  Nseg<-(I2-I0)
  Ic1<-Ic-I0
  oout<-list()
  if(Ic>I0){
    Y0<-Y[(I0+1):I2]
    T0<-T[(I0+1):I2]
    oout1<-LSmatrix(Y0,T0,NA)
    oout2<-LSmatrix(Y0,T0,Ic1)
    Fc<-(oout1$SSE-oout2$SSE)*(Nseg-3)/oout2$SSE
    prob<-pf(Fc,1,(Nseg-3))
    Pk0<-Pk.PMFT(Nseg)
    PFc<-Fc*Pk0[Ic1]
    oout$Fc<-Fc
    oout$PFc<-PFc
    oout$prob<-prob
  }
  else{
    oout$Fc<-0
    oout$PFc<-0
    oout$prob<-0
  }
  return(oout)
}

PMFxKxI0I2<-function(Y,T,I0,I2){
  Nmin2<-Nmin*2
  prob<-(-1)
  Ic<-I0
  Nseg<-(I2-I0)
  if(Nseg>=Nmin2){
    Y0<-Y[(I0+1):I2]
    T0<-T[(I0+1):I2]
    Pk0<-Pk.PMFT(Nseg)
    oout<-PMFT(Y0,T0,Pk0)
    Ic<-I0+oout$KPx
    prob<-pf(oout$Fx,1,(Nseg-3))
  }
  oout<-list()
  oout$Ic<-Ic
  oout$prob<-prob
  return(oout)
}

rmCycle<-function(idata){
  tdata<-cbind(idata[,2]*100+idata[,3],idata[,4])
  inds<-sort(unique(tdata[,1]))
  nx<-length(inds)
  mu<-rep(0,nx)
  for(i in 1:nx){
    mu[i]<-mean(tdata[tdata[,1]==inds[i],2],na.rm=T)
    tdata[tdata[,1]==inds[i],2]<-tdata[tdata[,1]==inds[i],2]-mu[i]
  }
  oout<-list()
  oout$EB<-mu
  oout$Base<-tdata[,2]
  return(oout)
}

LSmultiple<-function(Y,T,Ips){
  Nx<-length(Y)
  Ns<-length(Ips)-1
  X<-t(t(Y))
  D<-rep(1,Nx)
  D<-cbind(D,T-mean(T))
  if(Ns>=1){
    for(i in 1:Ns){
      tmp<-rep(0,Nx)
      tmp[(Ips[i]+1):Ips[i+1]]<-1
      D<-cbind(D,tmp)
    }
  }
  sig<-solve(t(D)%*%D)%*%t(D)%*%X
  fitted<-D%*%sig
  resi<-X-fitted
  SSE<-sum(resi^2)
  oout<-list()
  oout$SSE<-SSE
  oout$fitted<-as.vector(fitted)
  oout$resi<-as.vector(resi)
  oout$sig<-as.vector(sig)
  return(oout)
}

Rphi<-function(Y0,Ips,Ns){
# calculate auto-correlation of given data vector Y0 and breakpoints Ips
# output: cor -- autocorrelation; W -- prewhitenning vector of Y0
#  corl -- lower bound of cor; corh -- upper bound of cor
# if(Ns!=length(Ips)-1) stop("input data length error in Rphi")
  Y<-Y0
  N<-length(Y0)
  mu<-rep(0,Ns+1)
  for(i in 0:Ns){
    I1<- if(i==0) 1 else Ips[i]+1
    I2<-Ips[i+1]
    mu[i+1]<-mean(Y0[I1:I2])
    Y[I1:I2]<-Y0[I1:I2]-mu[i+1]
  }
  cor<-autocorlh(Y,IY0flg)
  W1<-Y
  W2<-Y
  W3<-Y
  W1[2:N]<-Y[2:N]-Y[1:(N-1)]*cor$cor
  W2[2:N]<-Y[2:N]-Y[1:(N-1)]*cor$corl
  W3[2:N]<-Y[2:N]-Y[1:(N-1)]*cor$corh
  W<-IY0flg*W1+(!IY0flg)*Y # if IY0flg==1 (continuous), W1; else Y
  WL<-IY0flg*W2+(!IY0flg)*Y # if IY0flg==1 (continuous), W1; else Y
  WU<-IY0flg*W3+(!IY0flg)*Y # if IY0flg==1 (continuous), W1; else Y
  for(i in 0:Ns){
    I1<- if(i==0) 1 else Ips[i]+1
    I2<-Ips[i+1]
    W[I1:I2]<-W[I1:I2]+mu[i+1]
    WL[I1:I2]<-WL[I1:I2]+mu[i+1]
    WU[I1:I2]<-WU[I1:I2]+mu[i+1]
  }
  oout<-list()
  oout$cor<-cor$cor
  oout$corl<-cor$corl
  oout$corh<-cor$corh
  oout$W<-W
  oout$WL<-WL
  oout$WU<-WU
  return(oout)
}

autocorlh<-function(Y,IY){
# calculate autocorrelation of given data vector, using given time vector to
# judge continuouse
  N<-length(Y)
  cnt<-sum(IY)
  m0<-mean(Y,na.rm=T)
  xsd0<-0
  xsd1<-0
  S1<-sum(((Y-m0)^2*IY)[1:(N-1)])
  S2<-sum(((Y-m0)*(c(Y[2:N],0)-m0)*IY)[1:(N-1)])
  cor<-S2/S1
# else stop("too few available data in autocor") 
  z975<-1.96
  z<-.5*log((1+cor)/(1-cor))
  df<-sum(IY[1:(N-1)])
  zl<-z-z975/sqrt(df-3)
  zh<-z+z975/sqrt(df-3)
  cl<-tanh(zl)
  ch<-tanh(zh)
  corl<-min(c(cl,ch))
  corh<-max(c(cl,ch))
  oout<-list()
  oout$cor<-cor
  oout$corl<-corl
  oout$corh<-corh
  return(oout)
}

getPFx95<-function(cor,N){
# if(cor<phi[1]|cor>phi[length(phi)]) stop("input series autocorrelation outbound!")
  if(cor<=phi[1])
    PTx95<-PFmax[N,1]
  else if(cor>=phi[length(phi)])
    PTx95<-PFmax[N,length(phi)]
  else{
    for(i in 1:(length(phi)-1))
      if(cor>phi[i]&cor<phi[i+1]) {
        Kr1<-i
        Kr2<-i+1
        cor1<-phi[i]
        cor2<-phi[i+1]
      }
    tmp1<-PFmax[N,Kr1]
    tmp2<-PFmax[N,Kr2]
    PTx95<-tmp1+(tmp2-tmp1)*(cor-cor1)/(cor2-cor1)
  }
  return(PTx95)
}

PMFxIseg<-function(Y0,Ti,Ips,Iseg){
  Ns<-length(Ips)-1
  N<-length(Y0)
  I0<- if(Iseg==1) 0 else Ips[Iseg-1]
  I3<-Ips[Iseg+1]
  Ic<-Ips[Iseg]
  Nseg<-I3-I0
  Ip0<-Ips[-Iseg]

  otmp<-LSmultipleRed(Y0,Ti,Ips)
  cor<-otmp$cor
  corl<-otmp$corl
  corh<-otmp$corh
  W<-otmp$W
  WL<-otmp$WL
  WU<-otmp$WU
  PFx95<-getPFx95(cor,Nseg)
  PFx95L<-getPFx95(corl,Nseg)
  PFx95U<-getPFx95(corh,Nseg)

  otmp<-LSmultiple(W,Ti,Ip0)
  SSE0.Iseg<-sum(otmp$resi[(I0+1):I3]^2)
  otmp<-LSmultiple(W,Ti,Ips)
  SSEf.Iseg<-sum(otmp$resi[(I0+1):I3]^2)
  Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
  if(Fx<0){
    Fx<-0
    PFx<-0
    prob<-0
  }
  else prob<-pf(Fx,1,(Nseg-3))

  otmp<-LSmultiple(WL,Ti,Ip0)
  SSE0.Iseg<-sum(otmp$resi[(I0+1):I3]^2)
  otmp<-LSmultiple(WL,Ti,Ips)
  SSEf.Iseg<-sum(otmp$resi[(I0+1):I3]^2)
  Fx1<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
  probL1<- if(Fx1<0) 0 else pf(Fx1,1,(Nseg-3))

  otmp<-LSmultiple(WU,Ti,Ip0)
  SSE0.Iseg<-sum(otmp$resi[(I0+1):I3]^2)
  otmp<-LSmultiple(WU,Ti,Ips)
  SSEf.Iseg<-sum(otmp$resi[(I0+1):I3]^2)
  Fx1<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
  probU1<- if(Fx1<0) 0 else pf(Fx1,1,(Nseg-3))

  probL<-min(c(probL1,probU1)); probU<-max(c(probL1,probU1))

  otmp<-LSmultiple(Y0,Ti,Ips)
# fitted<-otmp$fitted
  resi<-otmp$resi
# otmp<-Rphi(resi,Ips,Ns)
# W<-otmp$W+fitted
# cor<-otmp$cor
  SSEf.Iseg<-sum(resi[(I0+1):I3]^2)

  otmp<-LSmultiple(Y0,Ti,Ip0)
  SSE0.Iseg<-sum(otmp$resi[(I0+1):I3]^2)
  Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
  Pk0<-Pk.PMFT(Nseg)
  PFx<-Fx*Pk0[Ic-I0]

  oout<-list()
  oout$Fx<-Fx
  oout$PFx<-PFx
  oout$prob<-prob
  oout$probL<-probL
  oout$probU<-probU
  oout$PFx95<-PFx95
  oout$PFx95L<-PFx95L
  oout$PFx95U<-PFx95U
  return(oout)
}

LSmultipleRed<-function(Y0,Ti,Ips){
  Ns<-length(Ips)-1
  N<-length(Y0)
  otmp<-LSmultiple(Y0,Ti,Ips)
  sig<-otmp$sig
  beta<-otmp$sig[2]
  resi<-otmp$resi
  otmp<-autocorlh(resi,IY0flg)
  cor<-otmp$cor
  corl<-otmp$corl
  corh<-otmp$corh
  resi<-resi+beta*Ti
  W1<-resi/(1-cor)
  W2<-c(W1[1],(resi[2:N]-cor*resi[1:(N-1)])/(1-cor))
  W<-c(1,IY0flg[1:(N-1)])*W2+(!c(1,IY0flg[1:(N-1)]))*W1
  otmp<-LSmatrix(W,Ti,NA)
  beta<-otmp$sig[2]
  St0<-sum((Ti-mean(Ti))^2)
  df<-(N-2-Ns-Nt)
  sigmaE2<-otmp$SSE/df
  t.stat<-abs(beta)/sqrt(sigmaE2/St0)
  p.tr<-pt(t.stat,df)
  betaL<-beta-qt(.975,df)*sqrt(sigmaE2/St0)
  betaU<-beta+qt(.975,df)*sqrt(sigmaE2/St0)
  itmp<-Y0-beta*Ti
  mu<-rep(0,Ns+1)
  meanhat<-mu
  for(i in 0:Ns){
    I0<- if(i==0) 1 else Ips[i]+1
    I2<-Ips[i+1]
    mu[i+1]<-mean(itmp[I0:I2])
    meanhat[I0:I2]<-mu[i+1]+beta*Ti[I0:I2]
    resi[I0:I2]<-Y0[I0:I2]-meanhat[I0:I2]
  }
  W1<-resi
  W2<-c(resi[1],resi[2:N]-cor*resi[1:(N-1)])
  W3<-c(resi[1],resi[2:N]-corl*resi[1:(N-1)])
  W4<-c(resi[1],resi[2:N]-corh*resi[1:(N-1)])
  W<-c(1,IY0flg[1:(N-1)])*W2+(!c(1,IY0flg[1:(N-1)]))*W1
  WL<-c(1,IY0flg[1:(N-1)])*W3+(!c(1,IY0flg[1:(N-1)]))*W1
  WU<-c(1,IY0flg[1:(N-1)])*W4+(!c(1,IY0flg[1:(N-1)]))*W1
  for(i in 0:Ns){
    I0<- if(i==0) 1 else Ips[i]+1
    I2<-Ips[i+1]
    W[I0:I2]<-W[I0:I2]+mean(itmp[I0:I2])+beta*Ti[I0:I2]
    WL[I0:I2]<-WL[I0:I2]+mean(itmp[I0:I2])+beta*Ti[I0:I2]
    WU[I0:I2]<-WU[I0:I2]+mean(itmp[I0:I2])+beta*Ti[I0:I2]
  }
  oout<-list()
  oout$W<-W
  oout$WL<-WL
  oout$WU<-WU
  oout$sig<-sig
  oout$cor<-cor
  oout$corl<-corl
  oout$corh<-corh
  oout$resi<-resi
  oout$mu<-mu
  oout$meanhat<-meanhat
  oout$trend<-beta
  oout$betaL<-betaL
  oout$betaU<-betaU
  oout$p.tr<-p.tr
  return(oout)
}

LSmultiRedCycle<-function(Y0,Ti,Ips,Iseg.adj){
  N<-length(Y0)
  Ns<-length(Ips)-1
  Niter<-0
  tt<-TRUE
  EB1<-EB
  while(tt){
    tt<-FALSE
    Niter<-Niter+1
    EB0<-EB1
    otmp<-LSmultipleRed(Y0,Ti,Ips)
    trend<-otmp$trend
    betaL<-otmp$betaL
    betaU<-otmp$betaU
    resi<-otmp$resi
    cor<-otmp$cor
    corl<-otmp$corl
    corh<-otmp$corh
    p.tr<-otmp$p.tr
    meanhat<-otmp$meanhat
    mu<-otmp$mu
    W<-otmp$W
    WL<-otmp$WL
    WU<-otmp$WU

    if(Nt>1){
      itmp1<-cbind(EB0,Icy)
      itmp2<-cbind(1:N,Imd)
      colnames(itmp2)<-c("idx","Icy")
      itmp<-merge(itmp1,itmp2,by="Icy")
      EBfull<-itmp[order(itmp[,"idx"]),"EB0"]

      for(i in 1:(Ns+1)){
        I0<- if(i==1) 1 else Ips[i-1]+1
        I2<-Ips[i]
#       delta<-if(i==(Ns+1)) 0 else mu[i]-mu[Iseg.adj]
        delta<-mu[i]-mu[Iseg.adj]
        Y0[I0:I2]<-Y0[I0:I2]+EBfull[I0:I2]-delta
      }
    
      for(i in 1:Nt) EB1[i]<-mean(Y0[Imd==Icy[i]])
      VEB<-sqrt(var(EB1))
      if(is.na(VEB)) tt<-FALSE
      else{
        itmp1<-cbind(EB1,Icy)
        itmp2<-cbind(1:N,Imd)
        colnames(itmp2)<-c("idx","Icy")
        itmp<-merge(itmp1,itmp2,by="Icy")
        EBfull<-itmp[order(itmp[,"idx"]),"EB1"]

        for(i in 1:(Ns+1)){
          I0<- if(i==1) 1 else Ips[i-1]+1
          I2<-Ips[i]
#         delta<-if(i==(Ns+1)) 0 else mu[i]-mu[Iseg.adj]
          delta<-mu[i]-mu[Iseg.adj]
          Y0[I0:I2]<-Y0[I0:I2]-EBfull[I0:I2]+delta
        }

        DEBmx<-max(abs(EB1-EB0))
        if(DEBmx>VEB/1000&Niter<20) tt<-TRUE
      }
    }
  }
  oout<-list()
  oout$trend<-trend
  oout$betaL<-betaL
  oout$betaU<-betaU
  oout$EB<-EB1
  oout$mu<-mu
  oout$cor<-cor
  oout$corl<-corl
  oout$corh<-corh
  oout$W<-W
  oout$WL<-WL
  oout$WU<-WU
  oout$resi<-resi
  oout$Y0<-as.vector(Y0)
  oout$meanhat<-as.vector(meanhat)
  oout$p.tr<-p.tr
  return(oout)
}


QMadj.GaussianDLY.wRef<-function(Bseries,Rseries,MissingValueCode,InCs,output,Iadj=10000,Mq=20,Ny4a=5){
  if(Ny4a>0&Ny4a<=5) Ny4a<-5

  Read.wRef(Bseries,Rseries,MissingValueCode)
  N<-length(Y0); Nadj<-Nt*Ny4a
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

  ofileSout<-paste(output,"_QMadjDLYwRefstat.txt",sep="")
  file.create(ofileSout)
  cat(paste("The base data filename is:", Bseries,"; N=",N,"\n"),file=ofileSout)
  cat(paste("The ref data filename is:", Rseries,"\n"),file=ofileSout,append=T)

# if(Iadj==1|Iseg.longest==0) Iseg.adj<-Ns+1 else Iseg.adj<-Iseg.longest
  if(Iadj>(Ns+1)|Iseg.longest==0) Iseg.adj<-Ns+1
  else if(Iadj==0)Iseg.adj<-Iseg.longest
  else Iseg.adj<-Iadj

# estimate delta from Y0 (Base-Ref)
  otmp<-Rphi(Y0,Ips,Ns)
  cor<-otmp$cor
  muDif<-rep(0,Ns+1)
  for(i in 1:(Ns+1)){
    I0<- if(i==1) 1 else Ips[i-1]+1
    I2<- if(i>Ns) length(Y0) else Ips[i]
    muDif[i]<-mean(Y0[I0:I2])
  }

# transfer Ips(Base-Ref) to Ips(Base)
  Ips0<-Ips
  IY1<-bdata[,1]*10000+bdata[,2]*100+bdata[,3]
  IYM<-bdata[,2]*100+bdata[,3]
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
    adjBase[I0:I2]<-adjBase[I0:I2]+EB[I0:I2]-DeltaD # adjBase is base series adj to last seg
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
  dtmp<-LSmultipleRed(Aadj,Ti,Ipd)
# muD<-dtmp$mu[1]+EPBa-mean(bdata[,4])
  muD<-dtmp$mu[1]
  betaD<-dtmp$trend

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

  if(Ns>0){
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

  ofilePdf<-paste(output,"_QMadjDLYwRef.pdf",sep="")
  pdf(file=ofilePdf,onefile=T,paper='letter')
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  par(mfrow=c(3,1))
  par(mar=c(3,4,3,2)+.1,cex.main=.8,cex.lab=.8,cex.axis=.8,cex=.8)

  uyrs<-unique(floor(ori.bdata[,1]/10))*10
  labels<-NULL
  ats<-NULL
  for(i in 1:length(uyrs)){
    if(!is.na(match(uyrs[i],ori.bdata[,1]))){
      labels<-c(labels,uyrs[i])
      ats<-c(ats,match(uyrs[i],ori.bdata[,1]))
    }
  }

  pdata<-rep(NA,dim(ori.bdata)[1])
    pdata[owflg]<-Base+EB
  plot(1:dim(ori.bdata)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(Base+EB),max(Base+EB)),
       xaxt="n",col="black",lwd=0.5,
       main="a. Base series and Mean-adjusted  base series")
  axis(side=1,at=ats,labels=labels)
# pdata[owflg]<-otmp$meanhat+mean(EB1)
# lines(1:dim(ori.bdata)[1],pdata,col="red")
  pdata[owflg]<-meanhatD
  lines(1:dim(ori.bdata)[1],pdata,col="blue")
  pdata[owflg]<-adjB
  lines(1:dim(ori.bdata)[1],pdata,col="red")

  if(Ns>0) if(QMout$Mq>1){
    pdata[owflg]<-B
    plot(1:dim(ori.bdata)[1],pdata,type="l",xlab="",ylab="",
         ylim=c(min(B),max(B)),
         xaxt="n",col="black",lwd=0.5,
         main="b. QM-adjusted base series")
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

  if(Ns>0){
    odata<-matrix(NA,dim(ori.bdata)[1],8)
    odata[,1]<-ori.bdata[,1]
    odata[,2]<-ori.bdata[,2]
    odata[,3]<-ori.bdata[,3]
    odata[,4]<-ori.bdata[,4]
    if(QMout$Mq>1) odata[owflg,5]<-round(B,4)        # QMadjusted series
    odata[owflg,6]<-round(adjB,4)     # Mean-adjusted series
    if(QMout$Mq>1) odata[owflg,7]<-round(B-ori.bdata[owflg,4],4)
    odata[owflg,8]<-round(adjB-ori.bdata[owflg,4],4)
  }
  else odata<-cbind(ori.itable[,c(1,2,3,4,4,4)],0,0)
  ofileAout<-paste(output,"_QMadjDLYwRef.dat",sep="")
  write.table(odata,file=ofileAout,na=MissingValueCode,col.names=F,row.names=F)
  return(0)

}


PTKI0I2<-function(Y0,I0,I2){
# search new breakpoint of Y0[(I0+1):I2] using function PTK()
# output: Ic -- breakpoint, prob and PTx
  Y<-Y0[(I0+1):I2]
  N<-length(Y)
  oout<-list()
  oout$prob<-(-1)
  oout$Ic<-I0
  oout$PTx<-(-9999.9)
  if(N>=(Nmin*2)){
    Pk0<-Pk.PMT(N)
    otmp<-PTK(Y,Pk0)
    oout$Ic<-I0+otmp$KPx
    oout$prob<-pt(otmp$Tx,(N-2))
    oout$PTx<-otmp$PTx
  }
  return(oout)
}

PTK<-function(Y,Pk){
#  search input vector, return PTx.max and corresponding KPx
  PTx<-(-9999.9)
  KPx<-0
  N<-length(Y)
  for(k in Nmin:(N-Nmin)){
    EY1<-mean(Y[1:k])
    EY2<-mean(Y[(k+1):N])
    var<-sum(c((Y[1:k]-EY1)^2,(Y[(k+1):N]-EY2)^2))
    std<-sqrt(var/(N-2))
    Tk<-sqrt(k*(N-k)/N)*abs(EY1-EY2)/std
    PTk<-Tk*Pk[k]
    if(PTk>PTx){
      PTx<-PTk
      KPx<-k
      Tx<-Tk
    }
  }
  oout<-list()
  oout$PTx<-PTx
  oout$KPx<-KPx
  oout$Tx<-Tx
  return(oout)
}

Pk.PMT<-function(N){
# calculate penalty with given series length -- N
# output P, real vector of length N-1
  if(floor(N)!=N) stop("input para is not integer")
  Nle10<- if(N<=10) TRUE else FALSE
  Nlt50<- if(N<50) TRUE else FALSE
  Nle100<- if(N<=100) TRUE else FALSE
  Ngt10lt50<- if(N>10&N<50) TRUE else FALSE

  K<-seq(1,(N-1))
  A<-abs(1-2*K/N)
  B<-log(N)
  C<-log(B)
  D<-log(log(N+150))

  F <- if(Nle100) 1-A^((7*B-2*B*C)/10) else 1-A^(11*B*C/50)
  v <- if(Nle100) (15*C^0.5-11)/100 else (2*C^2+2*C-1)/100
  Po<-(11*C^(9/8)+195)*F^v/200

  K1<-sum(Po[1:floor(N/2)]<1)+1
  L<- if(Ngt10lt50) floor(K1/2)+2 else floor(K1/2)+1

  if(N<=10) delta<-rep(D^0.5*(Po[L+1]-Po[L]),(N-1))
  else if(N<=100) delta<-rep((D^(1/3)*(Po[L+1]-Po[L])+3/(10*N^(4/3))),(N-1))
  else{
    delta<-rep(0,(N-1))
    delta[1:L]<-(Po[L]-Po[1])*A[1:L]^(C^3)/(2*L-4)
    delta[(N-L):(N-1)]<-(Po[N-L]-Po[N-1])*A[(N-L):(N-1)]^(C^3)/(2*L-4)
  }
  P<-Po
  P[1:L]<-P[L]-delta[1:L]*(L-c(1:L))
  P[(N-L):(N-1)]<-P[N-L]-delta[(N-L):(N-1)]*(c((N-L):(N-1))-N+L)
  return(P)
}

PTKIc<-function(Y,Pk,Ic){
# calculate PTk and prob for given data vector and breakpoint Ic
  N<-length(Y)
  oout<-list()
  if(Ic>0){
    EY1<-mean(Y[1:Ic])
    EY2<-mean(Y[(Ic+1):N])
    var<-sum(c((Y[1:Ic]-EY1)^2,(Y[(Ic+1):N]-EY2)^2))
    std<-sqrt(var/(N-2))
    Tk<-sqrt(Ic*(N-Ic)/N)*abs(EY1-EY2)/std
    PTk<-Tk*Pk[Ic]
    oout$PTk<-PTk
    oout$prob<-pt(Tk,(N-2))
  }
  else{
    oout$PTk<-0
    oout$prob<-0
  }
  return(oout)
}

getPTx95<-function(cor,N){
# if(cor<phi[1]|cor>phi[length(phi)]) stop("input series autocorrelation outbound!")
  if(cor<=phi[1])
    PTx95<-PTmax[N,1]
  else if(cor>=phi[length(phi)])
    PTx95<-PTmax[N,length(phi)]
  else{
    for(i in 1:(length(phi)-1))
      if(cor>=phi[i]&cor<phi[i+1]) {
        Kr1<-i
        Kr2<-i+1
        cor1<-phi[i]
        cor2<-phi[i+1]
      }
    tmp1<-PTmax[N,Kr1]
    tmp2<-PTmax[N,Kr2]
    PTx95<-tmp1+(tmp2-tmp1)*(cor-cor1)/(cor2-cor1)
  }
  return(PTx95)
}


str40<-function(vari.name,vari.value){
  olen<-40
  b20<-"                    "
  if(!exists(vari.name,env=sys.parent())) 
    otmp<-paste(b20,"                  NA",sep="")
  else{
    if(nchar(vari.value)>olen) otmp<-paste("...",substr(vari.value,
       nchar(vari.value)-olen+4,nchar(vari.value)),sep="")
    else{
      otmp<-paste(b20,b20,sep="")
      substr(otmp,olen+1-nchar(vari.value),olen)<-vari.value
    }
  }
  return(otmp)
}


OnStepSize.wRef<-function(){
  getfile1<-function(){
    if(!exists("ifname")) ifname<-tclvalue(tkgetOpenFile())
    else ifname<-tclvalue(tkgetOpenFile(initialfile=ifname))
    if(!nchar(ifname)){
      tkinsert(txt,"end","No Base Data file selected in StepSize.wRef!\n\n")
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
      tkinsert(txt,"end","No Ref Data file selected in StepSize.wRef!\n\n")
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
  tkwm.title(tt,"StepSize.wRef")
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
      itmp<-StepSize.wRef(Bseries=ifname,Rseries=ifrname,InCs=iDIpsName,
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
        oact<-str40("ofbody",ofbody)
        oref<-str40("ifrname",ifrname)
        ocurdir<-str40("curdir",curdir)
        ooutdir<-str40("outdir",outdir)
        UIpsName1<-paste(ofname,"_fCs.txt",sep="")
        UIpsName<-paste(ofname,"_mCs.txt",sep="")
        file.copy(from=UIpsName1,to=UIpsName,overwrite=TRUE)
        tkgrid(tklabel(frameMiddle,text=oact,width=40),column=2,row=6,sticky="e")
        tkgrid(tklabel(frameMiddle,text=oref,width=40),column=2,row=7,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ocurdir,width=40),column=2,row=8,sticky="e")
        tkgrid(tklabel(frameMiddle,text=ooutdir,width=40),column=2,row=9,sticky="e")
        tkinsert(txt,"end",paste("data: ",ifname,"\noutput: ",ofname,"_*\n",sep=""))
        if(itmp<0) tkinsert(txt,"end",ErrorMSG)
        else{
          tkinsert(txt,"end","StepSize.wRef finished successfully...\n")
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

OnQuit<-function(){
  if(exists("curdir")) rm("curdir",envir=.GlobalEnv)
  if(exists("outdir")) rm("outdir",envir=.GlobalEnv)
  if(exists("ofbody")) rm("ofbody",envir=.GlobalEnv)
  if(exists("ofref")) rm("ofref",envir=.GlobalEnv)
  if(exists("ofname")) rm("ofname",envir=.GlobalEnv)
  if(exists("ifname")) rm("ifname",envir=.GlobalEnv)
  if(exists("ifrname")) rm("ifrname",envir=.GlobalEnv)
  if(exists("iIpsName")) rm("iIpsName",envir=.GlobalEnv)
  if(exists("iDIpsName")) rm("iDIpsName",envir=.GlobalEnv)
  if(exists("MissingStr")) rm("MissingStr",envir=.GlobalEnv)
  if(exists("PlevStr")) rm("PlevStr",envir=.GlobalEnv)
  if(exists("AdjStr")) rm("AdjStr",envir=.GlobalEnv)
  if(exists("Mq0Str")) rm("Mq0Str",envir=.GlobalEnv)
  if(exists("Ny4a")) rm("Ny4a",envir=.GlobalEnv)
  if(exists("xscr")) rm("xscr",envir=.GlobalEnv)
  if(exists("yscr")) rm("yscr",envir=.GlobalEnv)
  if(exists("txt")) rm("txt",envir=.GlobalEnv)
  if(exists("ErrorMSG")) rm("ErrorMSG",envir=.GlobalEnv)
  if(exists("textMissing")) rm("textMissing",envir=.GlobalEnv)
  tkdestroy(main)
  if(exists("main")) rm("main",envir=.GlobalEnv)
}

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
  tkgrid(tklabel(tt,text="Please enter the nominal conf. level p.lev value."),
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
      GuiErrorMSG<-paste(GuiErrorMSG,"p.lev must be one of these: 0.75,0.80,0.90,0.95,0.99,0.9999. Please re-enter\n")
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
      tkinsert(txt,"end",paste("The nominal level p.lev = ",PlevStr,"..\n",sep=" "))
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
