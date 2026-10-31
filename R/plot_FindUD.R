# plot_FindUD(oout, ofilePdf, EBfull, EEB, B, QMout, Ms, Mq, Ns, adj, Ips, Iseg.adj)
plot_FindUD <- function(oout, ofilePdf, EBfull, EEB, B, QMout, Ms, Mq, Ns, adj,
    Ips, Iseg.adj,
    otmp, ...)
{
    pdf(file=ofilePdf, onefile=TRUE, paper='letter')
    op <- par(no.readonly = TRUE) # the whole list of settable par's.
    par(mfrow=c(2,1))
    par(mar=c(3,4,3,2)+.1,cex.main=.8,cex.lab=.8,cex.axis=.8,cex=.8)

    uyrs<-unique(floor(ori.itable[,1]/10))*10
    labels<-NULL
    ats<-NULL
    for(i in 1:length(uyrs)){
        if(!is.na(match(uyrs[i],ori.itable[,1]))){
            labels<-c(labels,uyrs[i])
            ats<-c(ats,match(uyrs[i],ori.itable[,1]))
        }
    }

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
         ylim=c(min(otmp$Y0+EBfull,otmp$meanhat+EEB),
                max(otmp$Y0+EBfull,otmp$meanhat+EEB)),
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

    if(QMout$Mq>1){
        pdata[ooflg]<-B
        plot(1:dim(ori.itable)[1],pdata,type="l",xlab="",ylab="",
             ylim=c(min(c(adj,B)),max(c(adj,B))),
             xaxt="n",col="black",lwd=0.5,
             main="QM-adjusted base series")
        axis(side=1,at=ats,labels=labels)
    }

    # QM 订正量分布图
    if (Ns > 0 && QMout$Mq > 1) plot_qm_adjustments(QMout, Ips, Iseg.adj, Ns)

    par(op)
    dev.off()
}
