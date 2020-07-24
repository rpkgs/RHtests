plot_StepSize.ref <- function(output, Base, EB, EB1, B, 
                           oY0, omuDif, otmp, meanhatD, 
                           QMout, Mq, Ns, adj, adjB, Ips, Iseg.adj, ...)
{
  ofilePdf <- paste0(output, "_F.pdf")
  pdf(file=ofilePdf,onefile=T,paper='letter')
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  on.exit({ par(op); dev.off() })

  par(mfrow=c(4, 1))
  par(mar=c(3,4,3,2)+.1,mgp=c(1,.5,0),cex.main=.8,cex.lab=.8,cex.axis=.8,cex=.8)

  uyrs   <- unique(floor(ori.bdata$year / 10)) * 10 # decades
  ats    <- match(uyrs, ori.bdata$year) %>% rm_empty()
  labels <- ori.bdata$year[ats]

  IYori <- ori.bdata[,1]*10000+ori.bdata[,2]*100+ori.bdata[,3]
  rtmp  <- cbind(IY0,oY0,omuDif)
  stmp  <- merge(rtmp,t(t(IYori)),all.y=TRUE,by.x="IY0",by.y=1)
  pdata <- stmp[,2]

  plot(1:dim(ori.bdata)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(oY0),max(oY0)),
       xaxt="n",col="black",lwd=0.5,
       main="a. Base-minus-reference series")
  axis(side=1,at=ats,labels=labels)

  pdata<-stmp[,3]
  lines(1:dim(ori.bdata)[1],pdata,col="red")

  pdata[owflg]<-Base
  plot(1:dim(ori.bdata)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(Base,otmp$meanhat),max(Base,otmp$meanhat)),
       xaxt="n",col="black",lwd=0.5,
       main="b. De-seasonalized base series")
  axis(side=1,at=ats,labels=labels)
  pdata[owflg]<-otmp$meanhat
  lines(1:dim(ori.bdata)[1],pdata,col="red")

  pdata[owflg]<-Base+EB
  plot(1:dim(ori.bdata)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(Base+EB),max(Base+EB)),
       xaxt="n",col="black",lwd=0.5,
       main="c. Base series")
  axis(side=1,at=ats,labels=labels)
  pdata[owflg]<-otmp$meanhat+mean(EB1)
  lines(1:dim(ori.bdata)[1],pdata,col="red")
  pdata[owflg]<-meanhatD
  lines(1:dim(ori.bdata)[1],pdata,col="blue")

  xr<-dim(ori.bdata)[1]*.1; yrr<-mean((Base+EB)[1:xr])
  if(abs(yrr-max(Base+EB))<abs(yrr-min(Base+EB))) {
    yrr1<-min(Base+EB)+.2*(max(Base+EB)-min(Base+EB))
    yrr2<-min(Base+EB)+.05*(max(Base+EB)-min(Base+EB))
  } else {
    yrr1<-max(Base+EB)-.05*(max(Base+EB)-min(Base+EB))
    yrr2<-max(Base+EB)-.2*(max(Base+EB)-min(Base+EB))
  }
  text(2,yrr1,"Adjustments estimated from Base",adj=0,cex=.8); lines(dim(ori.bdata)[1]*c(.31,.33),rep(yrr1,2),col="red")
  text(2,yrr2,"Adjustments estimated from Diff",adj=0,cex=.8); lines(dim(ori.bdata)[1]*c(.31,.33),rep(yrr2,2),col="blue")

  pdata[owflg]<-adj
  plot(1:dim(ori.bdata)[1],pdata,type="l",xlab="",ylab="",
       ylim=c(min(c(adj,B)),max(c(adj,B))),
       xaxt="n",col="black",lwd=0.5,
       main="d. Mean-adjusted base series")
  axis(side=1,at=ats,labels=labels)
  pdata[owflg]<-adjB
  lines(1:dim(ori.bdata)[1],pdata,col="red")

  xr<-dim(ori.bdata)[1]*.1; yrr<-mean(adj[1:xr])
  if(abs(yrr-max(c(adj,B)))<abs(yrr-min(c(adj,B)))) {
    yrr1<-min(c(adj,B))+.2*(max(c(adj,B))-min(c(adj,B)))
    yrr2<-min(c(adj,B))+.05*(max(c(adj,B))-min(c(adj,B)))
  }
  else {
    yrr1<-max(c(adj,B))-.05*(max(c(adj,B))-min(c(adj,B)))
    yrr2<-max(c(adj,B))-.2*(max(c(adj,B))-min(c(adj,B)))
  }
  text(2,yrr1,"Adjustments estimated from Base",adj=0,cex=.8); lines(dim(ori.bdata)[1]*c(.31,.33),rep(yrr1,2),col="black")
  text(2,yrr2,"Adjustments estimated from Diff",adj=0,cex=.8); lines(dim(ori.bdata)[1]*c(.31,.33),rep(yrr2,2),col="red")

  if(Ns>0) if(QMout$Mq>1){
    pdata[owflg]<-B
    plot(1:dim(ori.bdata)[1],pdata,type="l",xlab="",ylab="",
         ylim=c(min(c(adj,B)),max(c(adj,B))),
         xaxt="n",col="black",lwd=0.5,
         main="e. QM-adjusted base series")
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
}
