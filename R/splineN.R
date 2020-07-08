
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
