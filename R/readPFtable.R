readPFtable<-function(Nx, pkth=4, ...){
  # ,ifile="PFmax31red_Nmin10_CVs"
  # read in PTmax table, assign PTmax table as global variable;
  # phi -- vector for cor catalog -- as global variable
  # itmp <- matrix(scan(ifile,skip=2,quiet=T),ncol=6,byrow=T)
  # itmp <- read.table("PFtable.csv", header = FALSE) %>% as.matrix()
  itmp <- PFtable
  # itmp<-itmp[,c(1,2,2,2,3,4,5,6)] # temp, will remove later
  Nmax  <- max(itmp[,1])
  Nmin  <- min(itmp[,1])
  ncol  <- dim(itmp)[2]
  nlevs <- unique(itmp[,1])
  PTmax <- matrix(0,max(Nmax,Nx),45)
  phi   <- c(-0.975,itmp[1:44,ncol])

  for(i in 0:(length(nlevs)-1)){
    for(j in 1:44){
      if(i==0){
        PTmax1 <- 0
        PTmax2 <- itmp[j,(pkth+1)]
        ind1   <- 1
        ind2   <- nlevs[1]
      } else {
        PTmax1 <- itmp[(i*44-44+j),(pkth+1)]
        PTmax2 <- itmp[(i*44+j),(pkth+1)]
        ind1   <- nlevs[i]
        ind2   <- nlevs[i+1]
      }
      PTmax[ind1:ind2,(j+1)] <- seq(PTmax1, PTmax2, length=(ind2-ind1+1))
    }
    PTmax[ind1:ind2,1] <- PTmax[ind1:ind2,2]-(phi[2]-phi[1])*
      (PTmax[ind1:ind2,3]-PTmax[ind1:ind2,2])/(phi[3]-phi[2])
  }
  if(Nx>Nmax)
    for(j in 1:45){
      delta<-(PTmax[nlevs[length(nlevs)],j]-PTmax[nlevs[length(nlevs)-1],j])/
        (nlevs[length(nlevs)]-nlevs[length(nlevs)-1])
      PTmax[Nmax:Nx,j]<-seq(from=PTmax[nlevs[length(nlevs)],j],
                            length=(Nx-Nmax+1),by=delta)
    }
  # PTmax<-PTmax[-1,]
  assign("phi", phi, envir=.GlobalEnv)
  assign("PFmax", PTmax, envir=.GlobalEnv)
  listk(phi, PFmax = PTmax)
}
