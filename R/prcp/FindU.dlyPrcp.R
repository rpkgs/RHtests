FindU.dlyPrcp<-function(InSeries, output, MissingValueCode="-999.99",
  GUI=FALSE,pthr=0.0, Mq=10, plev=0.95,Iadj=10000,Ny4a=0)
{
  # ErrorMSG<-NA
  # assign("ErrorMSG",ErrorMSG,envir=.GlobalEnv)
  Nmin<-10
  Ncat.min<-20
  # if(!p.lev%in%c(0.75,0.8,0.9,0.95,0.99,0.9999)){
  #   ErrorMSG<<-paste("FindU.dlyPrcp: input p.lev",p.lev,"error\n",
  #                    get("ErrorMSG",env=.GlobalEnv),"\n")
  #   if(!GUI) cat(ErrorMSG)
  #   return(-1)
  # }
  # pkth<-match(p.lev,c(0.75,0.8,0.9,0.95,0.99,0.9999))
  assign("Nmin",Nmin,envir=.GlobalEnv)
  itmp <- ReadDLY(InSeries,MissingValueCode,pthr)
  
  if(Ny4a>0&Ny4a<=5) Ny4a<-5
  Nall<-dim(ori.itable)[1]
  N<-length(Y0); Pfreq<-N/Nall; Nadj<-Ny4a*366*Pfreq
  critical_values <- readPFtable(N, plev)
  phi <- critical_values$phi
  PFmax <- critical_values$PFmax

  # smallestP<-min(Y0)
  smallestP<-min(ori.itable[ori.itable[,"data"]>0,"data"],na.rm=TRUE)
  smallestP<-max(c(smallestP,pthr))
  sumLogPi<-sum(log(Y0))

  ofileSout<-paste(output,"_Ustat.txt",sep="")
  file.create(ofileSout)
  fprintf("Input data filename: ; N= %s ; Nday= %s ; pthr= %s \n", N, Nall, pthr, file = ofileSout)
  fprintf("The nominal level of confidence (1-alpha)= %s \n", plev, file = ofileSout, append = TRUE)
  fprintf("The smallest non-zero dailyP is: %s \n", smallestP, file = ofileSout, append = TRUE)

  P       <- Y0
  aa      <- (-1)
  bb      <- 1
  lambdas <- seq(-1,1,.01)
  LHs     <- rep(NA,length(lambdas))
  Nss     <- rep(NA,length(lambdas))
  Ipss    <- matrix(NA,length(lambdas),100)
  Niter0  <- 0
  LH.max  <- (-9999999.)
  lr.max  <- (-1)
  for(dr in c(1,0.5,0.2,0.1)){
    lrs<-seq(aa,bb,by=dr)
    for(lr in 1:length(lrs)){
      lambda<-aa+(lr-1)*dr
      ind<-round((lambda+1)*100+1)
      if(abs(lambdas[ind]-lambda)>1e-8) {
        print(c(lambdas[ind],lambda,lambdas[ind]-lambda))
        stop(paste("ind=",ind,"lambda=",lambda,"error"))
      }
      if(is.na(LHs[ind])==FALSE){  # already done before
        LH<-LHs[ind]
        Ns<-Nss[ind]
        Ips<-Ipss[ind,1:(Ns+1)]
        Y1<-BCtrans(P,lambda)
      }
      else{
        Niter0<-Niter0+1
        Y1<-BCtrans(P,lambda)

        Ip0<-N
        Pk0<-Pk.PMFT(N)

        oout<-PMFT(Y1,Ti,Pk0, Nmin)
        I0<-0
        I2<-oout$KPx
        I4<-N
        I1<-PMFxKxI0I2(Y1,Ti,I0,I2, Nmin)$Ic
        I3<-PMFxKxI0I2(Y1,Ti,I2,I4, Nmin)$Ic
        I2<-PMFxKxI0I2(Y1,Ti,I1,I3, Nmin)$Ic

        Ns <- 1
        Ips <- c(I3, N)
        indices <- c(I1, I2, I3)
        candidates <- pmf_candidates(indices, Y1, Ti, N, IY0flg)
        PFx <- vapply(candidates, `[[`, numeric(1), "PFx")
        best <- sort(PFx, decreasing = TRUE, index.return = TRUE)$ix[1]
        PFx.mx <- PFx[best]
        prob.mx <- candidates[[best]]$prob
        Imx <- indices[best]
        cor.mx <- candidates[[best]]$corL
        PFx95L<-getPFx95(cor.mx,N, phi, PFmax)
        if(PFx.mx<PFx95L){
          Ns<-0
          Ips<-N
        }
        else{
          Ns<-1
          Ips<-c(Imx,N)

          tt<-TRUE
          Niter<-0
          while(tt){  # condition on is there any more Bps to insert in Ips?
            Niter<-Niter+1
            tt<-FALSE
            Ips0<-NULL
            for(i in 1:(Ns+1)){ # search for new break points
              I0<- if(i==1) 0 else Ips[i-1]
              I2<-Ips[i]
              otmp<-PMFxKxI0I2(Y1,Ti,I0,I2, Nmin)
              if(otmp$prob>0) Ips0<-sort(c(Ips0,otmp$Ic))
            }
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
                  otmp<-PMFxIseg(Y1,Ti,Ips1,ith, IY0flg, Nt, phi, PFmax)
                  if(otmp$PFx<otmp$PFx95L) Ips0[i]<-0
                  else if(otmp$PFx>PFx.mx){
                    PFx.mx<-otmp$PFx
                    Iseg.mx<-Ips0[i]
                    PFx95L.mx<-otmp$PFx95L
                  }
                }
                if(PFx.mx>=PFx95L.mx){ # Ic is significant
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

          # start to delete the least significant changepoint if it is insignificant
          tt<-TRUE
          while(tt){
            tt<-FALSE
            Iseg.mn<-0
            PFx.mn<-99999
            PFx95L.mn<-99999
            for(i in 1:Ns){
              otmp<-PMFxIseg(Y1,Ti,Ips,i, IY0flg, Nt, phi, PFmax)
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
        otmp<-LSmultipleRed(Y1,Ti,Ips, IY0flg, Nt)
        SSEf<-sum(otmp$resi**2)
        LH<--log(SSEf/(N-2-Ns))*N/2+(lambda-1)*sumLogPi
        #     ind<-floor((lambda+1)*100+1)
        #     if(lambdas[ind]!=lambda) stop(paste("ind=",ind,"lambda=",lambda,"error"))
        #     print(paste("##1##",round(lambda,2),Ns,Niter0,Ns,Ips[1],round(LH,2)))
        LHs[ind]<-LH
        Nss[ind]<-Ns
        if(Ns>0) Ipss[ind,1:(Ns+1)]<-Ips
      }
      if(LH>LH.max){
        LH.max<-LH
        Ips.max<-Ips
        Ns.max<-Ns
        lambda.max<-lambda
      }
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

  fprintf("Best lambda= %s ; Ns= %s ; likelihood= %s \n\n",
    round(lambda,2), Ns, round(LH.max,2), file = ofileSout, append = TRUE
  )
  Y1<-BCtrans(Y0,lambda)

  if(Ns>0) {
    Nsegs<-Ips-c(0,Ips[1:Ns])
    Iseg.longest<-sort(Nsegs,index=TRUE,decreasing=TRUE)$ix[1]
  }
  else Iseg.longest<-0

  if(Iadj>(Ns+1)|Iseg.longest==0) Iseg.adj<-Ns+1
  else if(Iadj==0)Iseg.adj<-Iseg.longest
  else Iseg.adj<-Iadj

  Ipd<-c(N)
  dtmp<-LSmultipleRed(Y1,Ti,Ipd, IY0flg, Nt)
  beta0<-dtmp$trend
  # cat(paste("Ignore changepoints -> trend0 =",round(beta0,6),
  #   "(p=", round(dtmp$p.tr,4), "); cor=", round(dtmp$cor,4),
  #   "(",round(dtmp$corl,4),",",round(dtmp$corh,4),")\n\n"),
  #   file=ofileSout,append=TRUE)

  otmp<-LSmultipleRed(Y1,Ti,Ips, IY0flg, Nt)
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
  #     file=ofileSout,append=TRUE)
  # cat(paste("#steps= ",Ns,"; trend=",round(otmp$trend,6),"unit/yr (p=",
  #           round(otmp$p.tr,4),"); lambda=",round(lambda,2),"; cor=",
  #           round(cor,4),"(",round(corL,4),",",round(corU,4),") p=",
  #           round(pcor,4),"\n"),
  #           file=ofileSout,append=TRUE)

  if(Ns>0) {
    fprintf("Step-sizes of the transfromed/original series:\n", file = ofileSout, append = TRUE)
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
      fprintf("%s %s transfered stepsize= %s original stepsize= %s \n",
        Ips[i], IY0[Ips[i]], round(stepsize,4), round(C[i],4), file = ofileSout, append = TRUE
      )
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
    fprintf("Nseg_shortest = %s ; Mq =  %s \n",
      QMout$Nseg.mn, QMout$Mq, file = ofileSout, append = TRUE
    )
    fprintf("\n Adjust to segment %s : from %s to %s \n",
      Iseg.adj, if(Iseg.adj==1) 1 else Ips[Iseg.adj-1]+1, Ips[Iseg.adj], file = ofileSout, append = TRUE
    )

    if(QMout$Mq>1){
      cat("#Fcat, DP (CDF and Differnces in category mean)\n",file=ofileSout,
          append=TRUE)
      write.table(round(QMout$osmean,4),file=ofileSout,append=TRUE,
                  row.names=FALSE,col.names=FALSE)
      for(i in 1:(Ns+1)){
        I1<-if(i==1) 1 else Ips[i-1]+1
        I2<-Ips[i]
        fprintf("Seg. %s: mean of QM-adjustments =%s\n",
          i, round(mean(QMout$W[I1:I2]),4), file = ofileSout, append = TRUE
        )
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
    #   odata[Imd,5]<-round(PA,4)
    #   odata[Imd,6]<-round(PdA,4)
    #   odata[Imd,7]<-round(PA-ori.itable[Imd,4],4)
    #   odata[Imd,8]<-round(PdA-ori.itable[Imd,4],4)
    odataP<-cbind(itable[,1],IY0,itable[,5],Pmeanhat,PA,PdA,PmeanhatA,
                  Y1,meanhat,Y,meanhatA,PA-itable[,5],PdA-itable[,5])
    otrend.ori<-getMtrendFdly(itable[,2:5])
    otrend.IBC<-getMtrendFdly(cbind(itable[,2:4],PdA))
    fprintf("Linear trend in the original monthly total P is trend0= %s mm/month\n",
      round(otrend.ori,5), file = ofileSout, append = TRUE
    )
    if(QMout$Mq>1){
      otrend.QM<-getMtrendFdly(cbind(itable[,2:4],PA))
      fprintf("Linear trend in QMadjusted monthly total P is trendQM= %s mm/month\n",
        round(otrend.QM,5), file = ofileSout, append = TRUE
      )
    }
    fprintf("Linear trend in mean adjusted(IBC) monthly total P is trend(IBC)= %s mm/month\n\n",
      round(otrend.IBC,5), file = ofileSout, append = TRUE
    )
  }
  else{
    #   odata<-cbind(ori.itable[,c(1,2,3,4,4,4)],0,0)
    odataP<-cbind(itable[,1],IY0,itable[,5],Pmeanhat,itable[,5],itable[,5],
                  Pmeanhat,Y1,meanhat,Y1,meanhat,0,0)
    odata<-cbind(1:nrow(ori.itable),ori.itable[,1]*10000+ori.itable[,2]*100+ori.itable[,3],
                 ori.itable[,4],ori.itable[,4],0,0,0,0,0,0,0,0)
    otrend.ori<-getMtrendFdly(itable[,2:5])
    fprintf("Linear trend in the original monthly total P is trend0= %s mm/month\n",
      round(otrend.ori,5), file = ofileSout, append = TRUE
    )
  }

  ofilePout<-paste(output,"_U.dat",sep="")
  # write.table(round(odataP,4),file=ofilePout,na=MissingValueCode,col.names=FALSE,row.names=FALSE)
  write.table(odata,file=ofilePout,na=MissingValueCode,col.names=FALSE,row.names=FALSE)

  ofilePdf<-paste(output,"_U.pdf",sep="")
  pdf(file=ofilePdf,onefile=TRUE)
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  par(mfrow=c(2,1))
  par(mar=c(2,4,3,1)+.1)
  par(cex=.8,cex.main=.8,cex.lab=.8,cex.axis=.8)
  if(Ns>0){
    #   par(mar=c(3,4,3,2)+.1)

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

      tflg<-is.na(p1data[I1:I2,3])==FALSE
      lines(c(I1:I2)[tflg],p1data[I1:I2,3][tflg],col="red",lwd=1.5)
    }
    for(Iseg in 1:Ns){ # plot Y1~meanhat
      #     I1<-max(c(1,Ips[Iseg]-180))
      #     Ic<-Ips[Iseg]
      #     I2<-min(c(N,Ips[Iseg]+180))
      I1<-max(c(1,(itable[Ips[Iseg],1]-180)))
      Ic<-itable[Ips[Iseg],1]
      I2<-min(c(Nall,(itable[Ips[Iseg],1]+180)))
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
  pdata=rep(NA,nrow(ori.itable))
  pdata[itable[,1]]<-P
  plot(1:nrow(ori.itable),pdata,type="l",col="black",ylab="prcp(mm)",
       main="Original dailyP>pthr series",xaxt="n")
  axis(side=1,at=ats,labels=yrs)
  # lines(1:nrow(ori.itable),odataP[,4],col="red")

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
      ymax<-max(osp[,2:3],na.rm=TRUE); ymin<-min(osp[,2:3],na.rm=TRUE)
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

  otmp<-LS_MultiBreaks(Y1,Ti,Ips)
  resi<-otmp$resi
  otmpW<-LS_MultiBreaks(W,Ti,Ips)
  resiW<-otmpW$resi
  otmpWL<-LS_MultiBreaks(WL,Ti,Ips)
  resiWL<-otmpWL$resi
  otmpWU<-LS_MultiBreaks(WU,Ti,Ips)
  resiWU<-otmpWU$resi

  ofileIout<-paste(output,"_1Cs.txt",sep="")
  ofileMout<-paste(output,"_mCs.txt",sep="")
  file.create(ofileIout)

  d_TP <- list()
  if(Ns==0) {
    fprintf("%s changepoints in Series \n", Ns, file = ofileIout)
    cat("PMF finds no Type-1 changepoints in the series!\n")
  }
  else{
    fprintf("%s changepoints in Series \n", Ns, file = ofileIout)
    
    for(i in 1:Ns){
      I1<- if(i==1) 1 else Ips[i-1]+1
      Ic<-Ips[i]
      I3<-Ips[i+1]
      Nseg<-I3-I1+1
      PFx95<-getPFx95(cor,Nseg, phi, PFmax)
      PFx95l<-getPFx95(corL,Nseg, phi, PFmax)
      PFx95h<-getPFx95(corU,Nseg, phi, PFmax)
      SSEf.Iseg<-sum(resi[I1:I3]^2)
      Ips1<-Ips[-i]
      otmp1<-LS_MultiBreaks(Y1,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      Pk1<-Pk.PMFT(Nseg)
      PFx<-Fx*Pk1[Ic-I1+1]

      SSEf.Iseg<-sum(resiW[I1:I3]^2)
      otmp1<-LS_MultiBreaks(W,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      Fx<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(Fx>0) prob<-pf(Fx,1,Nseg-3)
      else{
        Fx<-0
        prob<-0
        PFx<-0
      }

      SSEf.Iseg<-sum(resiWL[I1:I3]^2)
      otmp1<-LS_MultiBreaks(WL,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      FxL<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(FxL>0) probL<-pf(Fx,1,Nseg-3)
      else probL<-0

      SSEf.Iseg<-sum(resiWU[I1:I3]^2)
      otmp1<-LS_MultiBreaks(WU,Ti,Ips1)
      SSE0.Iseg<-sum(otmp1$resi[I1:I3]^2)
      FxU<-(SSE0.Iseg-SSEf.Iseg)*(Nseg-3)/SSEf.Iseg
      if(FxU>0) probU<-pf(Fx,1,Nseg-3)
      else probU<-0

      #     else if(Id==1) { # type-1 changepoints
      if(PFx<PFx95l) Idc<-"No  "
      else if(PFx>=PFx95l&PFx<PFx95h) Idc<-"?   "
      else if(PFx>=PFx95h) Idc<-"Yes "
      #     }

      fprintf("%s %s%s (%s-%s)%s%s (%s-%s)\n",
        sprintf("%1.0f",1), sprintf("%-4.4s",Idc), sprintf("%10.0f", IY0[Ic]),
        sprintf("%10.4f",probL), sprintf("%10.4f",probU), sprintf("%6.3f",plev),
        sprintf("%10.4f",PFx), sprintf("%10.4f",PFx95l), sprintf("%10.4f",PFx95h),
        file = ofileIout, append = TRUE
      )
      fprintf("PMF : c=%s; (Time %s); Type= 1; p=%s(%s-%s); PFmax=%s; CV95=%s(%s-%s); Nseg=%s\n",
        sprintf("%4.0f",Ic), sprintf("%10.0f",IY0[Ic]), sprintf("%10.4f",prob),
        sprintf("%10.4f",probL), sprintf("%10.4f",probU), sprintf("%10.4f",PFx),
        sprintf("%10.4f",PFx95), sprintf("%10.4f",PFx95l), sprintf("%10.4f",PFx95h),
        sprintf("%4.0f",Nseg), file = ofileSout, append = TRUE
      )
      d_TP[[i]] <- data.table(kind = 1, Idc, date = IY0[Ic],
            probL, probU, plev, PFx, PFx95l, PFx95h)
    }

    if (!is_empty(d_TP)) {
        d_TP %<>% do.call(rbind, .)
        d_TP[, 4:9] %<>% lapply(round, digits = 4)
        fwrite(d_TP, ofileIout, append = TRUE, col.names = TRUE)
    }
  }
  
  if(GUI)
    return(0)
  else{
    file.copy(from=ofileIout,to=ofileMout,overwrite=TRUE)
    cat("FindU.dlyPrcp finished successfully...\n")
    # odata
    list(data = odata, TP = d_TP)
  }
}
