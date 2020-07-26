
plot_stepsize <- function(oout, ofilePdf, EBfull, EEB, B, QMout, Ms, Mq, Ns, adj,
    Ips, Iseg.adj,
    otmp, ...)
{
    pdf(file=ofilePdf,onefile=T,paper='letter')
    op <- par(no.readonly = TRUE) # the whole list of settable par's.
    par(mfrow=c(2,1),cex.main=.8,cex.lab=.8,cex.axis=.8,cex=.8)

    uyrs   <- unique(floor(ori.itable[,1]/10))*10
    labels <- NULL
    ats    <- NULL
    for(i in 1:length(uyrs)){
      if(!is.na(match(uyrs[i],ori.itable[,1]))){
        labels <- c(labels,uyrs[i])
        ats    <- c(ats,match(uyrs[i],ori.itable[,1]))
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
          } else {
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
