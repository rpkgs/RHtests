
ReadDLY <- function(InSeries, MissingValueCode = "-999.99", pthr = 0.0){
    itable <- if (is.character(InSeries)) {
        sep = ifelse(is.csv(InSeries), ",", " ")
        fread(InSeries, sep = sep, na.strings = MissingValueCode)
    } else InSeries
    itable %<>% set_colnames(c("year","month","day","data")) %>% as.data.table()

    # keep input base data as ori.itable
    ori.itable <- itable
    Nall       <- dim(itable)[1]
    ind        <- c(1:Nall)
    otable     <- cbind(I = ind, itable)[!is.na(data) & data > pthr]
    Nday       <- otable[, I]
    IY0        <- otable[, year*10000+month*100+day]
    IY0flg     <- rep(0,length(IY0))
    Y0         <- otable[, data]
    Iyr        <- otable[, year]
    Imd        <- otable[, month*100+day]
    Ti         <- Nday/365.25

    for(i in 1:(length(IY0)-1)){
        if(Nday[i+1]-Nday[i]==1) IY0flg[i]<-1
    }
    if(sum(IY0flg)<10){  # too few available data for autocorlh
        ErrorMSG<<-paste("Too many missing values in ", idata, "to estimate autocorrelation\n")
        return(-1)
    }

    itable <- otable
    assign("ori.itable", as.data.frame(ori.itable), envir = .GlobalEnv)
    assign("itable", as.data.frame(itable), envir = .GlobalEnv)

    assign("Ti",Ti,envir=.GlobalEnv) # Time index for LS fitting
    assign("Y0",Y0,envir=.GlobalEnv) # Data series for Base
    assign("IY0",IY0,envir=.GlobalEnv) # Cycle index for Base
    assign("Imd",Imd,envir=.GlobalEnv) # Cycle index for Base
    assign("IY0flg",IY0flg,envir=.GlobalEnv) # continuous flag for Base
    d = data.table(IY0, Imd, IY0flg, Ti, Y0)
    d
    # return(0)
}
