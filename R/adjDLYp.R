
adjDLYp<-function(P,Ips,Mq,Iseg.adj,Ptr.mx,R,Ncat.min,Nadj){
  Ns<-length(Ips)-1
  N<-length(P)
  Nseg.mn<-N
  for(i in 1:(Ns+1)){
    I1<-if(i==1) 1 else Ips[i-1]+1
    I2<-Ips[i]
    Nseg<-I2-I1+1
    if(Nseg<Nseg.mn) Nseg.mn<-Nseg
  }
  if(Mq<=0) Mq<-min(floor(Nseg.mn/Ncat.min),100)
  else Mq<-min(floor(Nseg.mn/Ncat.min),Mq)
  if(Mq>100) Mq<-100
  if(Mq<=0) Mq<-1
  Fd<-.5/Mq
  Fcat<-matrix(NA,(Ns+1),(Mq+2))
  F<-matrix(NA,(Ns+1),N)
  EPb<-matrix(NA,(Ns+1),(Mq+2))
  EPa<-matrix(NA,(Ns+1),(Mq+2))
  for(i in 1:(Ns+1)) Fcat[i,]<-seq(0,by=1/Mq,length=(Mq+2))
  for(i in 1:Ns){
    I1<-if(i==1) 0 else Ips[i-1]
    I2<-Ips[i]
    if(Nadj>0) I1<-max(c(I1,I2-Nadj))
    Nseg<-I2-I1
    Y<-P[(I1+1):I2]
    iindx<-sort(Y,index=T)$ix
    irank<-sort(iindx,index=T)$ix
    F[i,1:Nseg]<-irank/Nseg
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
      #     F[i,1:Nseg]<-irank/Nseg
      EPa[i,]<-0
      for(k in 2:(Mq+1)){
        Mp1<-floor(Nseg*Fcat[i,(k-1)])
        Mp2<-floor(Nseg*Fcat[i,k])
        EPa[i,k]<-mean(Y[iindx[(Mp1+1):Mp2]])
      }
      EPa[i,(Mq+2)]<-EPa[i,(Mq+1)]
    }
  }
  # if(Mq==1) EPb[,1]<-EPb[,2]
  if(Nadj==0) if(Ns>1) EPa[1:(Ns-1),]<-EPb[2:Ns,]
  Adj<-matrix(0,(Ns+1),(Mq+2))
  for(k in 1:(Mq+2)){
    if(Iseg.adj>1) for(i in 1:(Iseg.adj-1)) Adj[i,k]<-sum(EPa[i:(Iseg.adj-1),k])-sum(EPb[i:(Iseg.adj-1),k])
    if(Iseg.adj<=Ns) for(i in (Iseg.adj+1):(Ns+1)) Adj[i,k]<-sum(EPb[Iseg.adj:(i-1),k])-sum(EPa[Iseg.adj:(i-1),k])
  }

  osmean<-c(1:(Mq+2)) # for plot purpose
  for(Iseg in c(1:(Ns+1)))
    osmean<-cbind(osmean,Fcat[Iseg,]-Fd,Adj[Iseg,])
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
    if(i==Iseg.adj) PA[I1:I2]<-P[I1:I2]-Ptr.mx+R[I1:I2] else{
      dx<-Fcat[i,]-Fd
      #     fdx<-EP[Iseg.adj,]-EP[i,]
      fdx<-Adj[i,]
      if(Mq==1) fdx[1]<-fdx[2]
      fdx2<-splineN(dx,fdx,2E30,2E30)
      for(j in I1:I2) W[j]<-splintN(dx,fdx,fdx2,F[i,(j-I1+1)])
      PA[I1:I2]<-P[I1:I2]+W[I1:I2]-Ptr.mx+R[I1:I2]
    }

    Rs<-F[i,1:Nseg]
    ors<-sort(Rs,index=T)$ix
    osp<-rbind(osp,cbind(I1:I2,Rs[ors],W[I1:I2][ors]))
  }

  oout<-list()
  oout$PA<-PA
  oout$W<-W
  oout$Nseg.mn<-Nseg.mn
  oout$Mq<-Mq
  oout$osmean<-osmean
  oout$osp<-osp
  return(oout)
}
