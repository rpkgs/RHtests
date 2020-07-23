#' QMadj
#'
#' @return
#'
#' Example1_Ref1_DLY.dat_QMadjDLY_data.dat, which contains the original and the
#' adjusted daily data in its 4 th and 5 th column respectively.
#'
#' (column 6, 7, 8 are the mean-adjusted daily series, the QM adjustments,
#' and the mean adjustments, respectively;
#'
#' columns 4 + 7 = column 5;
#' columns 4 + 8 = column 6)
#'
#' @export
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

#' @rdname QMadjGaussian
#' @export
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

        if(i==Iseg.adj) {
          BA[I1:I2]<-B[I1:I2]
          osp<-rbind(osp,cbind(I1:I2,0,0))
        }
        else{
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
            ttmp=cbind(Rs[ors],DB[ors])
            for(qth in 1:Mq)
              osmean[qth+1,i*2+1]=mean(ttmp[ttmp[,1]>=(qth-1)/Mq&ttmp[,1]<=qth/Mq,2])
            osmean[1,i*2+1]=osmean[2,i*2+1]
            osmean[Mq+2,i*2+1]=osmean[Mq+1,i*2+1]

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
