#  read-in data from base_data_file and ref_data_file, put MissingValueCode as NA
#  set several global variables as output:
#  Nt     --  total numbers of seasonal factors, say, monthly=12, daily=365
#  Icy    --  catelog of seasonal factors, say, monthly 1:12
#  Y0     --  data vector of base - ref, take common period for both non-missing
#  IY0    --  cor-responding seasonal vector of Y0
#  IY0flg --  flg(integer, 0 or 1) for IY0, 1 -- continouse, 0 -- not continouse
#             for autocorrelation calculation
#  bdata  --  matrix of non-missing base data, 4 columns, yyyy,mm,dd,data
#  ori.bdata  -- original base data matrix, also 4 columns, same as bdata
Read.wRef <- function(Bseries, Rseries, MissingValueCode="-999.99"){
    if (is.null(Bseries) || is.null(Rseries)) return()

    sep = ifelse(is.csv(Bseries), ",", "")

    var_colnames <- c("year", "month", "day", "data")
    # var_colnames <- c("id1", "id2", "id3", "data")
    itable <- if (is.character(Bseries)) {
        read.table(Bseries, sep = sep, header = FALSE,
            na.strings = MissingValueCode, colClasses = "real")
    } else Bseries
    
    rtable <- if (is.character(Rseries)) {
        read.table(Rseries, sep = sep, header = FALSE,
            na.strings = MissingValueCode, colClasses = "real")
    } else Rseries

    # check input data (both base and ref), if column!=4, return error
    itable %<>% set_colnames(var_colnames)
    rtable %<>% set_colnames(var_colnames)

    # ooflg      <- itable[, which(isnot_0229 & !is.na(data))]

    # # get rid of Feb 29th data
    # Nall <- n <- nrow(itable)
    # itable   <- itable[day != 29, ]
    # iyrbegin <- itable$year[1]
    # imdbegin <- itable$month[1]*100 + itable$day[1]
    # iyrend   <- itable$year[n]
    # imdend   <- itable$month[n]*100 + itable$day[n]

    # browser()
    # keep input base data as ori.itable
    ori.itable <- itable
    isnot_0229 <- itable[, !(month == 2 & day == 29)]
    owflg <- itable[, !is.na(data) & !(month == 2 & day == 29)]
    # get rid of Feb 29th data
    itable <- itable[isnot_0229, ]
    rtable <- rtable[isnot_0229, ]
    # check input data (both base and ref), no jump with begin and end
    Icy <- sort(unique(itable[,2]*100+itable[,3]))
    Nt  <- length(Icy)
    # construct YYYYMMDD for base series
    imdbegin <- itable[1,2]*100+itable[1,3] # begin MMDD for base series
    iyrbegin <- itable[1,1] # begin year for base series
    Nx1<-dim(itable)[1]
    imdend<-itable[Nx1,2]*100+itable[Nx1,3] # end MMDD for base series
    iyrend<-itable[Nx1,1] # end year for base series
    Ind1<-iyrbegin*10000+Icy[Icy>=imdbegin] # first year
    if(iyrend>(iyrbegin+1)) for(i in (iyrbegin+1):(iyrend-1))
        Ind1<-c(Ind1,i*10000+Icy)
    Ind1<-c(Ind1,iyrend*10000+Icy[Icy<=imdend])
    YMD.base<-itable[,1]*10000+itable[,2]*100+itable[,3]
    for(i in 1:length(Ind1)) if(Ind1[i]!=YMD.base[i]|is.na(YMD.base[i]))
        stop(paste("input base series not continuous at:",Ind1[i],YMD.base[i]))

    # construct YYYYMMDD for ref series
    imdbegin<-rtable[1,2]*100+rtable[1,3] # begin MMDD for ref series
    iyrbegin<-rtable[1,1] # begin year for base series
    Nx2<-dim(rtable)[1]
    imdend<-rtable[Nx2,2]*100+rtable[Nx2,3] # end MMDD for ref series
    iyrend<-rtable[Nx2,1] # end year for ref series
    Ind2<-iyrbegin*10000+Icy[Icy>=imdbegin] # first year
    if(iyrend>(iyrbegin+1)) for(i in (iyrbegin+1):(iyrend-1))
        Ind2<-c(Ind2,i*10000+Icy)
    Ind2<-c(Ind2,iyrend*10000+Icy[Icy<=imdend])
    YMD.ref<-rtable[,1]*10000+rtable[,2]*100+rtable[,3]
    for(i in 1:length(Ind2)) if(Ind2[i]!=YMD.ref[i]|is.na(YMD.ref[i])) {
        ErrorMSG<-paste("input ref series not continuous at:",Ind2[i],YMD.base[i],
                        "\n",get("ErrorMSG",env=.GlobalEnv))
        return(-1)
    }

    # take non-missing data only
    itable<-itable[is.na(itable[,4])==F,]
    # itable.nm<-itable
    rtable<-rtable[is.na(rtable[,4])==F,]

    Nx1        <- dim(itable)[1]
    Nx2        <- dim(rtable)[1]
    icol.base  <- itable[,1]*10000+itable[,2]*100+itable[,3]
    icol.ref   <- rtable[,1]*10000+rtable[,2]*100+rtable[,3]

    itable.nmb <- merge(cbind(icol.base,itable), cbind(icol.ref,rtable),by=1,all.x=T,all.y=F,sort=F)[,c(1,2,3,4,5,9)]
    itable.nmb <- itable.nmb[order(itable.nmb[,1]),]
    itable.nmb <- itable.nmb[,-1]
    colnames(itable.nmb) <- c(colnames(itable),'data.ref')

    ind.base<-cbind(icol.base,seq(1,Nx1))
    ind.ref<-cbind(icol.ref,seq(1,Nx2))
    ind.base<-ind.base[is.na(itable[,4])==F,]
    ind.ref<-ind.ref[is.na(rtable[,4])==F,]
    colnames(ind.base)<-c("IY0","ind")
    colnames(ind.ref)<-c("IY0","ind")
    cind<-merge(ind.base,ind.ref,by.x="IY0",by.y="IY0",
                suffixes=c(".base",".ref"))
    IY0<-cind[,"IY0"]
    IY0flg<-rep(0,length(IY0))
    # construct flag vector for autocor calculation
    Iyr<-floor(IY0/10000)
    Imd<-IY0-Iyr*10000
    Ti<-IY0
    for(i in 1:length(IY0)){
        ith<-match(Imd[i],Icy)
        Ti[i]<-(Iyr[i]-iyrbegin)*Nt+ith
    }
    IyrB<-floor(icol.base/10000)
    ImdB<-icol.base-IyrB*10000
    TiB<-rep(0,Nx1)
    for(i in 1:Nx1){
        ith<-match(ImdB[i],Icy)
        TiB[i]<-(IyrB[i]-iyrbegin)*Nt+ith
    }
    for(i in 1:(length(IY0)-1)){
        if(Ti[i+1]-Ti[i]==1) IY0flg[i]<-1
    }
    IYBflg<-rep(0,length(TiB))
    for(i in 1:(length(TiB)-1))
        if(TiB[i+1]-TiB[i]==1) IYBflg[i]<-1
    ind.base<-cind[,"ind.base"]
    ind.ref<-cind[,"ind.ref"]
    # check data qualification
    itmp<-cbind(itable[,2]*100+itable[,3],rep(NA,dim(itable)[1]))
    itmp[ind.base,2]<-itable[ind.base,4]
    idenind<-unique(itmp[,1])
    for(i in 1:Nt){
        if(sum(is.na(itmp[itmp[,1]==Icy[i]])==F)<10){
            ErrorMSG<<-paste("input data too few at:",Icy[i],
                             "\n",get("ErrorMSG",env=.GlobalEnv))
            return(-1)
        }
    }
    itmp1<-0
    for(i in 1:(dim(itmp)[1]-1))
        if(is.na(itmp[i,2])==F&is.na(itmp[(i+1),2])==F) itmp1<-itmp1+1
    if(itmp1<10) {
        ErrorMSG<<-paste("input data too few for autocorrelation calculating",
                         "\n",get("ErrorMSG",env=.GlobalEnv))
        return(-1)
    }
    # finish checking
    Y0   <- itable[ind.base,4]-rtable[ind.ref,4]
    rtmp <- itable[ind.base,]
    otmp <- rmCycle(rtmp)
    EBb  <- otmp$EB
    rtmp <- rtable[ind.ref,]
    otmp <- rmCycle(rtmp)
    EBr  <- otmp$EB
    itmp <- itable[ind.base,2]*100+itable[ind.base,3]
    for(i in 1:length(Y0)){
        indd<-itmp[i] # mmdd for Y0[i]
        indf<-NULL
        for(j in 1:Nt) if(Icy[j]==indd) indf<-j
        Y0[i]<-Y0[i]+EBr[indf]-EBb[indf]
    }

    browser()
    # c("id", "date", "value", "is_continue")
    listk(Ti, IY0, Y0, IY0flg) %>% str()
    listk(bdata = itable.nmb, ori.bdata = ori.itable, TiB, IYBflg, owflg, Icy, Nt) %>% str()

    assign("Ti", Ti, envir = .GlobalEnv) # Time index for LS fitting
    assign("Y0", Y0, envir = .GlobalEnv) # Data series for Base-Ref
    assign("IY0", IY0, envir = .GlobalEnv) # Cycle index for Base-Ref
    assign("IY0flg", IY0flg, envir = .GlobalEnv) # continuous flag for Base-Ref

    assign("bdata", itable.nmb, envir = .GlobalEnv) # non-missing table for base data
    assign("ori.bdata", ori.itable, envir = .GlobalEnv) # original base data
    assign("TiB", TiB, envir = .GlobalEnv)
    assign("IYBflg", IYBflg, envir = .GlobalEnv) # continuous flag for Base-Ref
    assign("owflg", owflg, envir = .GlobalEnv)

    assign("Icy", Icy, envir = .GlobalEnv) # Cycle index
    assign("Nt", Nt, envir = .GlobalEnv)   # Cycle length
}
