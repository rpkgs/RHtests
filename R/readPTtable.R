readPTtable<-function(Nx, plev)
{
  if (!plev %in% c(0.75, 0.8, 0.9, 0.95, 0.99, 0.9999)) {
    stop(paste0("invalid : input plev = ", plev, " !\n"))
  }
  pkth <- match(plev, c(0.75, 0.8, 0.9, 0.95, 0.99, 0.9999))

  # read in PTmax table, assign PTmax table as global variable;
  # phi -- vector for cor catalog -- as global variable
  # itmp <- read.table("PTtable.csv", header = FALSE) %>% as.matrix()
  itmp   <- PTtable
  # itmp <- matrix(, ncol=8,byrow=T)
  Nmax   <- max(itmp[,1])
  Nmin   <- min(itmp[,1])
  nlevs  <- unique(itmp[,1])
  PTmax  <- matrix(0,max(Nmax,Nx),45)
  phi    <- c(-0.975,itmp[1:44,8])
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
      PTmax[ind1:ind2,(j+1)]<-seq(PTmax1,PTmax2,length=(ind2-ind1+1))
    }
    PTmax[ind1:ind2,1]<-PTmax[ind1:ind2,2]-(phi[2]-phi[1])*
      (PTmax[ind1:ind2,3]-PTmax[ind1:ind2,2])/(phi[3]-phi[2])
  }
  if(Nx>Nmax)
    for(j in 1:45){
      delta <- (PTmax[nlevs[length(nlevs)],j]-PTmax[nlevs[length(nlevs)-1],j])/
        (nlevs[length(nlevs)]-nlevs[length(nlevs)-1])
      PTmax[Nmax:Nx,j]<-seq(from=PTmax[nlevs[length(nlevs)],j],
                            length=(Nx-Nmax+1),by=delta)
    }
  # PTmax<-PTmax[-1,]
  assign("phi", phi, envir=.GlobalEnv)
  assign("PTmax", PTmax, envir=.GlobalEnv)
  listk(phi, PTmax)
}
