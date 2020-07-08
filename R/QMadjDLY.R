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
