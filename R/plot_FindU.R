
# global variales:
# - `ori.itable`
#
plot_FindU <- function(oout, ofilePdf, EBfull, EEB, B, QMout, Ms, Mq, Ns, adj,
    Ips, Iseg.adj)
{
    pdf(file=ofilePdf, onefile=TRUE, paper='letter')
    op <- par(no.readonly = TRUE) # the whole list of settable par's.
    par(mfrow=c(2,1))
    par(mar=c(3,4,3,2)+.1,cex.main=.8,cex.lab=.8,cex.axis=.8,cex=.8)
    on.exit({ par(op); dev.off() })

    uyrs   <- unique(floor(ori.itable$year/10))*10 # decades
    ats    <- match(uyrs, ori.itable$year) %>% rm_empty()
    labels <- ori.itable$year[ats]

    I_valid <- ooflg
    plot2(I_valid, oout$Y0, type="l",xlab="",ylab="",
         ylim= range2(oout$Y0, oout$meanhat),
         xaxt="n",col="black",lwd=0.5,
         main="Base anomaly series and regression fit",
         at=ats,labels=labels)
    lines(I_valid, oout$meanhat, col = "red")

    plot2(I_valid, oout$Y0+EBfull, type="l",xlab="",ylab="",
         ylim=range2(oout$Y0+EBfull,oout$meanhat+EEB),
         xaxt="n",col="black",lwd=0.5,
         main="Base series and regression fit",
         at=ats,labels=labels)
    lines(I_valid, oout$meanhat+EEB, col="red")

    plot2(I_valid, adj, type="l",xlab="",ylab="",
         ylim=c(min(c(adj,B)),max(c(adj,B))),
         xaxt="n",col="black",lwd=0.5,
         main="Mean-adjusted base series",
         at=ats, labels=labels)

    if(Mq>1){
        plot2(I_valid, B, type="l",xlab="",ylab="",
             ylim=c(min(c(adj,B)),max(c(adj,B))),
             xaxt="n",col="black",lwd=0.5,
             main="QM-adjusted base series",
             at=ats,labels=labels)
    }

    # test plot
    if(Ns>0 & Mq>1){
        par(mar=c(4,5,3,2)+.1, cex=.8, mfrow=c(2,2), mgp=c(1.2,.5,0))
        col = 0
        np     <- 0
        osp    <- QMout$osp
        osmean <- QMout$osmean
        for(i in 1:(Ns+1)){
            Fd   <- .5/QMout$Mq
            I1   <- if(i==1) 1 else Ips[i-1]+1
            I2   <- Ips[i]
            ymax <- max(osp[,2:3],na.rm=T); ymin<-min(osp[,2:3],na.rm=T)
            if(i!=Iseg.adj){
                np <- np+1
                if(col==0) {
                    col <- 2
                    plot(osp[I1:I2,2], osp[I1:I2,3], xlim=c(0,1), ylim=c(ymin,ymax),
                         type="l",lwd=1,col=col,xlab="Cumulative Frequency",
                         ylab="QM Adjustment")
                    title(cex.main=0.9,main=paste("distribution of QM adjustments with Mq=",QMout$Mq),line=.5)
                    icol <- 2*np
                    for(j in 1:QMout$Mq){
                        lines(c(osmean[(j+1),icol]-Fd,osmean[(j+1),icol]+Fd),
                              c(rep(osmean[(j+1),(icol+1)],2)),col=col,lty=2,lwd=.5)
                        if(j>=1 & j<QMout$Mq)
                            lines(rep(osmean[(j+1),icol]+Fd,2),
                                c(osmean[(j+1),(icol+1)],osmean[(j+2),(icol+1)]),col=col,lty=2,lwd=0.5)
                    }
                } else {
                    col<-col+1
                    lines(osp[I1:I2,2],osp[I1:I2,3],lwd=1,col=col)
                    icol<-2*np
                    for(j in 1:QMout$Mq){
                        lines(c(osmean[(j+1),icol]-Fd,osmean[(j+1),icol]+Fd),
                              c(rep(osmean[(j+1),(icol+1)],2)),col=col,lty=2,lwd=.5)
                        if(j>=1 & j<QMout$Mq)
                            lines(rep(osmean[(j+1),icol]+Fd,2),
                                c(osmean[(j+1),(icol+1)],osmean[(j+2),(icol+1)]),col=col,lty=2,lwd=0.5)
                    }
                }
                text(.15,ymax-np*(ymax-ymin)/(Ns*3),paste("Seg.",i))
                lines(c(.25,.30),rep(ymax-np*(ymax-ymin)/(Ns*3),2),lwd=2,col=col)
            }
            else np<-np+1
        }
    }
}
