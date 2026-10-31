
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

    # QM 订正量分布图
    if (Ns > 0 && Mq > 1) plot_qm_adjustments(QMout, Ips, Iseg.adj, Ns)
}
