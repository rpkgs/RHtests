
FindUD.dlyPrcp<-function(InSeries, InCs, output,
  MissingValueCode="-999.99", GUI=FALSE, pthr=0.0,Mq=10,plev=0.95,Iadj=10000,Ny4a=0)
{
  # ErrorMSG<-NA
  # assign("ErrorMSG",ErrorMSG,envir=.GlobalEnv)
  Nmin<-10
  if(Ny4a>0&Ny4a<=5) Ny4a<-5
  Ncat.min<-20
  # if(!p.lev%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
  #   ErrorMSG<<-paste("FindU.dlyPrcp: input p.lev",p.lev,"error\n",
  #                    get("ErrorMSG",env=.GlobalEnv),"\n")
  #   if(!GUI) cat(ErrorMSG)
  #   return(-1)
  # }
  # plev<-p.lev
  # pkth<-match(p.lev,c(0.75,0.8,0.9,0.95,0.99,0.9999))
  assign("Nmin",Nmin,envir=.GlobalEnv)
  itmp <- ReadDLY(InSeries, MissingValueCode, pthr)
  # if(itmp==(-1)){
  #   ErrorMSG<<-paste("FindUD.dlyprcp: Error in read data from",InSeries,"\n",
  #                    get("ErrorMSG",env=.GlobalEnv),"\n")
  #   if(!GUI) cat(ErrorMSG)
  #   return(-1)
  # }

  Nall<-dim(ori.itable)[1]
  N<-length(Y0); Pfreq<-N/Nall; Nadj<-Ny4a*366*Pfreq
  readPFtable(N, plev)

  # itmp<-readLines(InCs)
  TP = if (is.character(InCs)) fread(InCs) else InCs
  if (is.Date(TP$date)) TP$date %<>% date2num()
  Ns <- nrow(TP) # number of changing points
  if (is.null(Ns)) Ns <- 0

  if(Ns == 0) {
    Ns.ini  <- 0
    Ips.ini <- N
    Ids.ini <- 0
  } else {
    Ns.ini <- Ns
    Ips.ini<-c(rep(0,Ns.ini),N)
    Ids.ini<-rep(0,Ns.ini)
    for(i in 1:Ns.ini){ # using YYYYMMDD as index, searching for the largest
      # date less or equal to given YYYYMMDD
      # ymdtmp<-as.numeric(substr(itmp[i+1],7,16))
      ymdtmp <- TP$date[i]
      it<-match(ymdtmp,IY0)
      if (!is.na(it)) {
        Ips.ini[i] <- it
      } else {
        Ips.ini[i] <- max(c(1:N)[IY0 <= ymdtmp])
      }
      # Ids.ini[i]<-as.numeric(substr(itmp[i+1],1,1))
    }
    Ids = TP$kind

    if(sum(is.na(Ips.ini))>0|!identical(Ips.ini,sort(Ips.ini))){
      ErrorMSG<<-paste("FindUD.dlyprcp: Ips read in from ",InCs,"error:")
      for(i in 1:Ns.ini)
        ErrorMSG<<-paste(get("ErrorMSG",env=.GlobalEnv),Ips[i])
      ErrorMSG<<-paste(get("ErrorMSG",env=.GlobalEnv),"\n\n")
      if(!GUI) cat(ErrorMSG)
      return(-1)
    }
  }

  smallestP <- min(ori.itable[ori.itable[,"data"]>0,"data"],na.rm=T)
  smallestP <- max(c(smallestP,pthr))
  sumLogPi  <- sum(log(Y0))

  ofileSout<-paste(output,"_UDstat.txt",sep="")
  file.create(ofileSout)
  cat(paste("Input data filename:","; N=",N,"; Nday=",Nall,
            "; pthr=",pthr,"\n"), file=ofileSout)
  cat(paste("The nominal level of confidence (1-alpha)=",plev,"\n"),
      file=ofileSout,append=T)
  cat(paste("The smallest non-zero dailyP is:",smallestP,"\n"),
      file=ofileSout,append=T)

  P<-Y0
  aa<-(-1)
  bb<-1
  lambdas<-seq(-1,1,.01)
  LHs<-rep(NA,length(lambdas))
  Nss<-rep(NA,length(lambdas))
  Ipss<-matrix(NA,length(lambdas),500)
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
        Ns<-Nss[ind]
        Ips<-Ipss[ind,1:(Ns+1)]
        Y1<-BCtrans(P,lambda)
      }
      else{
        tt<-TRUE
        Ns<-Ns.ini
        Ips<-Ips.ini
        Niter0<-Niter0+1
        Y1<-BCtrans(P,lambda)

        if(Ns==0){
          Ip0<-N
          Pk0<-Pk.PMFT(N)

          oout<-PMFT(Y1,Ti,Pk0)
          I0<-0
          I2<-oout$KPx
          I4<-N
          I1<-PMFxKxI0I2(Y1,Ti,I0,I2)$Ic
          I3<-PMFxKxI0I2(Y1,Ti,I2,I4)$Ic
          I2<-PMFxKxI0I2(Y1,Ti,I1,I3)$Ic

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
            corU1<-0
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
            corU2<-0
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
            corU3<-0
          }

          tmp<-sort(c(PFx1,PFx2,PFx3),decreasing=T,index.return=T)
          PFx.mx<-tmp$x[1]
          prob.mx<-c(prob1,prob2,prob3)[tmp$ix[1]]
          Imx<-c(I1,I2,I3)[tmp$ix[1]]
          cor.mx<-c(corL1,corL2,corL3)[tmp$ix[1]]
          PFx95L<-getPFx95(cor.mx,N)
          if(prob.mx<plev){
            Ns<-0
            tt<-FALSE
            Ids<-0
          }
          else{
            Ns<-1
            Ips<-c(Imx,N)
            Ids<-0
          }
        }
        else{  # there are input changepoints
          Ns<-Ns.ini
          Ips<-Ips.ini
          Ids<-Ids.ini
        }
        Ips.i<-Ips
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
          # finished search for possible new changepoints, start estimate the p-value
          # of the Ips0 series and find the most significant changepoint Iseg.mx
          tt1<-TRUE
          while(tt1){
            PFx.mx<-(-9999)
            Iseg.mx<-0
            prob.mx<-(-99999)
            probL.mx<-(-99999)
            if(length(Ips0)==0) tt1<-FALSE
            else{
              for(i in 1:length(Ips0)){
                Ips1<-sort(c(Ips,Ips0[i]))
                ith<-match(Ips0[i],Ips1)
                otmp<-PMFxIseg(Y1,Ti,Ips1,ith)
                probL<-min(c(otmp$probL,otmp$probU,otmp$prob))
                probU<-max(c(otmp$probL,otmp$probU,otmp$prob))
                #	      if(otmp$PFx<otmp$PFx95L) Ips0[i]<-0
                if(probU<plev) Ips0[i]<-0
                else if(otmp$PFx>PFx.mx){
                  prob.mx<-otmp$prob
                  probL.mx<-probL
                  Iseg.mx<-Ips0[i]
                  PFx.mx<-otmp$PFx
                }
              }
              if(probL.mx>=plev){ # Ic is significant
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
        # finish finding any possible new changepoints

        Ids0<-rep(NA,length(Ips))
        for(i in 1:length(Ips)){
          if(Ips[i]%in%Ips.i) Ids0[i]<-Ids[Ips.i==Ips[i]]
          else Ids0[i]<-0
        }
        Ids<-Ids0
        # start to delete the least significant changepoint if it is insignificant
        while(tt){
          tt<-FALSE
          Iseg.mn<-0
          probL.mn<-99999
          for(i in 1:Ns){
            if(Ids[i]==0){ # check those un-documented
              otmp<-PMFxIseg(Y1,Ti,Ips,i)
              probL<-min(otmp$probL,otmp$probU)
              if(probL<probL.mn){
                Iseg.mn<-i
                probL.mn<-probL
              }
            }
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
      if(length(Ips)==1) if(is.na(Ips)) Ips<-N

      otmp<-LSmultipleRed(Y1,Ti,Ips)
      SSEf<-sum(otmp$resi**2)
      LH<--log(SSEf/(N-2-Ns))*N/2+(lambda-1)*sumLogPi
      LHs[ind]<-LH
      Nss[ind]<-Ns
      if(Ns>0) Ipss[ind,1:(Ns+1)]<-Ips

      if(LH>LH.max){
        LH.max<-LH
        Ips.max<-Ips
        Ids.max<-Ids
        Ns.max<-Ns
        lambda.max<-lambda
      }
      #   cat(paste("ind=",ind,"lambda=",lambda,"Ns=",Ns,"likelihood=",round(LH,1),
      #             "lambda.max=",lambda.max,"\n"))
    }
    # end of detection
    if(dr>0.1){
      aa<-max(c(lambda.max-dr,-1))
      bb<-min(c(lambda.max+dr,1))
      lrs<-seq(aa,bb,by=dr)
    }
  }

  # final output
  lambda<-lambda.max
  Ns<-Ns.max
  Ips<-Ips.max
  # Ids<-Ids.max
  if(Ns>0){
    Ids0<-rep(NA,length(Ips))
    for(i in 1:length(Ips)){
      if(Ips[i]%in%Ips.ini&(i!=Ns)) Ids0[i]<-Ids.ini[which(Ips.ini==Ips[i])]
      else Ids0[i]<-0
    }
    Ids<-Ids0
  }
  else Ids<-0


  cat(paste("Best lambda=",round(lambda,2),"; Ns=",Ns,"; likelihood=",
            round(LH.max,2),"\n\n"),file=ofileSout,append=T)
  Y1<-BCtrans(Y0,lambda)

  if(Ns>0) {
    Nsegs<-Ips-c(0,Ips[1:Ns])
    Iseg.longest<-sort(Nsegs,index=T,decreasing=T)$ix[1]
  }
  else Iseg.longest<-0

  if(Iadj>(Ns+1)|Iseg.longest==0) Iseg.adj<-Ns+1
  else if(Iadj==0)Iseg.adj<-Iseg.longest
  else Iseg.adj<-Iadj

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
  Pmeanhat<-IVBCtrans(meanhat,lambda)
  meanhatA<-rep(NA,N)

  # cat("Common trend TPR fit to the transformed Base series:\n",
  #     file=ofileSout,append=T)
  # cat(paste("#steps= ",Ns,"; trend=",round(otmp$trend,6),"unit/yr (p=",
  #           round(otmp$p.tr,4),"); lambda=",round(lambda,2),"; cor=",
  #           round(cor,4),"(",round(corL,4),",",round(corU,4),") p=",
  #           round(pcor,4),"\n"),
  #           file=ofileSout,append=T)

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
    otrend.IBC<-getMtrendFdly(cbind(itable[,2:4],PdA))
    cat(paste("Linear trend in the original monthly total P is trend0=",
              round(otrend.ori,5),"mm/month\n"),file=ofileSout,append=T)
    if(QMout$Mq>1){
      otrend.QM<-getMtrendFdly(cbind(itable[,2:4],PA))
      cat(paste("Linear trend in QMadjusted monthly total P is trendQM=",
                round(otrend.QM,5),"mm/month\n"),file=ofileSout,append=T)
    }
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
  # ofileDout<-paste(output,"_adj_UD.dat",sep="")
  # write.table(odata,file=ofileDout,na=MissingValueCode,col.names=F,row.names=F)

  ofilePout<-paste(output,"_UD.dat",sep="")
  # write.table(round(odataP,4),file=ofilePout,na=MissingValueCode,col.names=F,row.names=F)
  write.table(odata,file=ofilePout,na=MissingValueCode,col.names=F,row.names=F)

  ofilePdf<-paste(output,"_UD.pdf",sep="")
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
      lines(c(I1:I2),odata[I1:I2,9],col="red",lwd=1.5)
    }
  }

  yrs<-unique(itable[,2])
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

  ofileIout<-paste(output,"_pCs.txt",sep="")
  ofileMout<-paste(output,"_mCs.txt",sep="")
  file.create(ofileIout)

  d_TP = list()
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

      if(Ids[i]==0) { # type-0 changepoints
        if(probU<plev) Idc<-"No  "
        else if(probL<plev&probU>=plev) Idc<-"?   "
        else if(probL>=plev) Idc<-"YifD"
        if(PFx>=PFx95h) Idc<-"Yes "
      }
      else if(Ids[i]==1){ # type-1 changepoints
        if(PFx<PFx95l) Idc<-"No  "
        else if(PFx>=PFx95l&PFx<PFx95h) Idc<-"?   "
        else if(PFx>=PFx95h) Idc<-"Yes "
      }

      cat(paste(sprintf("%1.0f",Ids[i])," ",
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
       d_TP[[i]] <- data.table(kind = Ids[i], Idc, date = IY0[Ic],
                probL, probU, plev, PFx, PFx95l, PFx95h) 
    }
    d_TP %<>% do.call(rbind, .)
    d_TP[, 4:9] %<>% lapply(round, digits = 4)
    fwrite(d_TP, ofileIout, append = TRUE, col.names = TRUE)
  }
  if(GUI)
    return(0)
  else{
    file.copy(from=ofileIout,to=ofileMout,overwrite=TRUE)
    cat("FindUD.dlyprcp finished successfully...\n")
    odata = NULL
    list(fit = odata, TP = d_TP)
  }
}
