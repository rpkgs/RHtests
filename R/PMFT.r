
Pk.PMFT<-function(N){
  # if(floor(N)!=N) stop("input data error in Pk")
  Nlt40<- if(N<40) TRUE else FALSE
  Nle500<- if(N<=500) TRUE else FALSE
  
  K<-seq(1,(N-1))
  Kmin<- if(floor((N-1)/2)==(N-1)/2) c(1:floor((N-1)/2),floor((N-1)/2):1)
         else c(1:(floor((N-1)/2)+1),floor((N-1)/2):1)
  W <- floor(if(Nle500) 11*N/50 else 21*N/100)
  A <- abs(1-2*K/N)
  B <- log(N)
  C <- log(B)
  D <- log(log(N+150))
  Q <- c(abs(1-50*K[1:floor(N/2)]/(N*11)),
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

PMFT<-function(Y,T,Pk0){
  N     <- length(Y)
  PFx   <- (-99999.)
  Fx    <- (-99999.)
  KPx   <- 0
  T_anorm = T - mean(T)
  oout1 <- LSmatrix(Y, T_anorm, NA)

  for(i in Nmin:(N-Nmin)){
    oout2 <- LSmatrix(Y, T_anorm, i, only.SSE = TRUE)
    Fc    <- (oout1$SSE-oout2$SSE)*(N-3)/oout2$SSE
    PFc   <- Fc*Pk0[i]
    if(PFc > PFx){
      PFx <- PFc
      Fx  <- Fc
      KPx <- i
    }
  }

  oout<-list()
  oout$PFx <- PFx
  oout$KPx <- KPx
  oout$Fx  <- Fx
  return(oout)
}

PMFxKc<-function(Y,T,I0,I2,Ic){
  Nseg<-(I2-I0)
  Ic1<-Ic-I0
  oout<-list()
  if(Ic>I0){
    Y0      <- Y[(I0+1):I2]
    T0      <- T[(I0+1):I2]

    T_anorm <- T0 - mean(T0)
    oout1   <- LSmatrix(Y0, T_anorm, NA)
    oout2   <- LSmatrix(Y0, T_anorm, Ic1)
    Fc      <- (oout1$SSE-oout2$SSE)*(Nseg-3)/oout2$SSE
    prob    <- pf(Fc,1,(Nseg-3))
    Pk0     <- Pk.PMFT(Nseg)
    PFc     <- Fc*Pk0[Ic1]
    oout$Fc   <- Fc
    oout$PFc  <- PFc
    oout$prob <- prob
  }
  else{
    oout$Fc   <- 0
    oout$PFc  <- 0
    oout$prob <- 0
  }
  return(oout)
}

PMFxKxI0I2<-function(Y,T,I0,I2){
  Nmin2 <- Nmin*2
  prob  <- (-1)
  Ic    <- I0
  Nseg  <- (I2-I0)
  if(Nseg>=Nmin2){
    Y0   <- Y[(I0+1):I2]
    T0   <- T[(I0+1):I2]
    Pk0  <- Pk.PMFT(Nseg)
    oout <- PMFT(Y0, T0, Pk0)
    Ic   <- I0+oout$KPx
    prob <- pf(oout$Fx,1,(Nseg-3))
  }
  oout<-list()
  oout$Ic   <- Ic
  oout$prob <- prob
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
  z975 <- 1.96
  z    <- .5*log((1+cor)/(1-cor))
  df   <- sum(IY[1:(N-1)])
  zl   <- z-z975/sqrt(df-3)
  zh   <- z+z975/sqrt(df-3)
  cl   <- tanh(zl)
  ch   <- tanh(zh)
  corl <- min(c(cl,ch))
  corh <- max(c(cl,ch))
  oout<-list()
  oout$cor<-cor
  oout$corl<-corl
  oout$corh<-corh
  return(oout)
}

getPFx95<-function(cor,N){
  n = length(phi)
# if(cor<phi[1]|cor>phi[length(phi)]) stop("input series autocorrelation outbound!")
  if(cor<=phi[1])
    PTx95<-PFmax[N,1]
  else if(cor>=phi[n])
    PTx95<-PFmax[N, n]
  else{
    Kr1 <- which(cor > phi[-n] & cor < phi[-1])
    Kr2 <- Kr1 + 1
    # for(i in 1:(length(phi)-1))
    #   if(cor>phi[i]&cor<phi[i+1]) {
    #     Kr1<-i
    #     Kr2<-i+1
    #   }
    cor1 <- phi[Kr1]
    cor2 <- phi[Kr2]
    tmp1 <- PFmax[N,Kr1]
    tmp2 <- PFmax[N,Kr2]
    PTx95 <- tmp1+(tmp2-tmp1)*(cor-cor1)/(cor2-cor1)
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

  otmp   <- LSmultipleRed(Y0,Ti,Ips)
  cor    <- otmp$cor
  corl   <- otmp$corl
  corh   <- otmp$corh
  W      <- otmp$W
  WL     <- otmp$WL
  WU     <- otmp$WU
  PFx95  <- getPFx95(cor,Nseg)
  PFx95L <- getPFx95(corl,Nseg)
  PFx95U <- getPFx95(corh,Nseg)

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
