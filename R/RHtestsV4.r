
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
