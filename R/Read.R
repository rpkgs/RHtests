#' read input of RHtestsV4
#'
#' Remove NA values and flag `y[t]` as 0 when it has not corresponding `y[t+1]` values.
#'
#' @inheritParams FindU
#'
#' @return
#' If failed, NULL will be returned. If success, a data.table with the following
#' columns will be returned:
#' - `IY0`: numeric, `year*1000 + month*100 + day`, e.g. '19980100'
#' - `IY0flg`: 0 or 1, is `y[t]` has `y[t+1]` for autocorrelation?
#' - `Imd`: numeric, `month*100 + day`, e.g. '100'
#' - `Ti`: index or original Y
#' - `Y0`: original Y after eliminating NA values
#'
#' `NT`=1, 12, and 365 for annual, monthly, and daily series, respectively.
#' @examples
#' infile = system.file("extdata/Example1.dat", package = "RHtests")
#' Read(infile, "-999.99")
#' @export
Read <- function(InSeries, MissingValueCode = "-999.99", plev = 0.95){    

    itable <- if (is.character(InSeries)) {
        sep = ifelse(is.csv(InSeries), ",", " ")
        fread(InSeries, sep = sep, na.strings = MissingValueCode)
    } else InSeries
    itable %<>% set_colnames(c("year","month","day","data")) %>% as.data.table()

    # keep input base data as ori.itable
    ori.itable <- itable
    isnot_0229 <- itable[, !(month == 2 & day == 29)]
    ooflg      <- itable[, which(isnot_0229 & !is.na(data))]

    # get rid of Feb 29th data
    
    itable   <- itable[isnot_0229, ]
    iyrbegin <- itable$year[1]
    imdbegin <- itable$month[1]*100 + itable$day[1]
    n <- nrow(itable)
    iyrend   <- itable$year[n]
    imdend   <- itable$month[n]*100 + itable$day[n]

    # check input data (both base and ref), no jump with begin and end
    Icy <- itable[, sort(unique(month * 100 + day))]
    # Icy  <- sort(unique(ori.itable[ooflg, 2]*100+ori.itable[ooflg,3]))
    Ind2 <- iyrbegin*10000 + Icy[Icy>=imdbegin] # first year

    if (iyrend > (iyrbegin + 1)) {
        for (i in (iyrbegin + 1):(iyrend - 1)) {
            Ind2 <- c(Ind2, i * 10000 + Icy)
        }
    }
    Ind2 <- c(Ind2,iyrend*10000+Icy[Icy<=imdend])
    Nt   <- length(Icy)

    ind  <- ori.itable[ooflg, year*10000+month*100+day]
    if(sum(!(ind%in%Ind2))>0|all.equal(ind,sort(ind))!=T){
        ErrorMSG <<- paste("Input data series not continuous\n")
        return(-1)
    } 
    # for(i in 1:length(Ind2))
    #   if(Ind2[i]!=ind[i]) {
    #     ErrorMSG<<-paste("Input data series not continuous at:",Ind2[i],ind[i],"\n")
    #     return(-1)
    #   }
    # IY0<-ind[is.na(itable[,4])==F]
    IY0    <- ind
    IY0flg <- rep(0, length(IY0))
    Y0     <- itable[!is.na(data), data]
    Iyr    <- floor(IY0/10000)
    Imd    <- IY0-Iyr*10000
    Ti     <- IY0
    for(i in 1:length(IY0)){
        ith<-match(Imd[i],Icy)
        Ti[i]<-(Iyr[i]-iyrbegin)*Nt+ith
    }

    for(i in 1:(length(IY0)-1)){
        if(Ti[i+1]-Ti[i]==1) IY0flg[i]<-1
    }
    if(sum(IY0flg)<10){  # too few available data for autocorlh
        ErrorMSG<<-paste("Too many missing values in ", "to estimate autocorrelation\n")
        return(-1)
    }

    N <- length(Y0)
    readPFtable(N, plev)

    itable <- itable[!is.na(data), ]
    assign("ori.itable", as.data.frame(ori.itable), envir = .GlobalEnv)
    assign("itable", as.data.frame(itable), envir = .GlobalEnv)

    assign("IY0", IY0, env = .GlobalEnv)
    assign("ooflg", ooflg, envir = .GlobalEnv) # numeric index, other than boolean
    assign("Ti", Ti, envir = .GlobalEnv)   # Time index for LS fitting
    assign("Y0", Y0, envir = .GlobalEnv)   # Data series for Base
    assign("IY0", IY0, envir = .GlobalEnv) # Cycle index for Base
    assign("Imd", Imd, envir = .GlobalEnv) # Cycle index for Base
    assign("IY0flg", IY0flg, envir = .GlobalEnv) # continuous flag for Base
    assign("Icy", Icy, envir = .GlobalEnv) # Cycle index
    assign("Nt", Nt, envir = .GlobalEnv)   # Cycle length
    d = data.table(IY0, Imd, IY0flg, Ti, Y0)
    d
    # return(0)
}

is.csv<-function(names){
    nlen<-nchar(names)
    if(substr(names,nlen-2,nlen)=="csv"|substr(names,nlen-2,nlen)=="CSV") return(T)
    else return(F)
}


ReadDLY.g<-function(InSeries,MissingValueCode){
    if(!file.exists(InSeries)) {
        ErrorMSG<<-paste("Input datafile",InSeries,"does not exist!\n")
        return(-1)
    }
    if(is.csv(InSeries)){
        itmp<-try(read.table(InSeries,sep=",",header=F,na.strings=MissingValueCode,
                             colClasses=rep("real",4)),silent=T)
        if(inherits(itmp,"try-error")){
            ErrorMSG<<-geterrmessage()
            return(-1)
        } else itable<-itmp
    }
    else{
        itmp<-try(read.table(InSeries,sep="",header=F,na.strings=MissingValueCode,
                             colClasses=rep("real",4)),silent=T)
        if(inherits(itmp,"try-error")){
            ErrorMSG<<-geterrmessage()
            return(-1)
        }
        else itable<-itmp
    }
    if(ncol(itable)!=4){
        ErrorMSG<<-paste(InSeries,"has",ncol(itable),"columns. The number of columns should be 4\n")
        return(-1)
    }
    colnames(itable)<-c("id1","id2","id3","data")

    iyrbegin<-itable[1,1]
    imdbegin<-itable[1,2]*100+itable[1,3]
    iyrend <- itable[dim(itable)[1],1]
    imdend <- itable[dim(itable)[1],2]*100+itable[dim(itable)[1],3]
    # keep input base data as ori.itable
    ori.itable<-itable
    # check input data (both base and ref), no jump with begin and end
    Icy<-sort(unique(itable[,2]*100+itable[,3]))
    Ind2<-iyrbegin*10000+Icy[Icy>=imdbegin] # first year
    # if(iyrend>(iyrbegin+1)) for(i in (iyrbegin+1):(iyrend-1))
    #   Ind2<-c(Ind2,i*10000+Icy)
    # Ind2<-c(Ind2,iyrend*10000+Icy[Icy<=imdend])
    Nt<-length(Icy)
    Nall<-dim(itable)[1]
    ind<-itable[,1]*10000+itable[,2]*100+itable[,3]
    # for(i in 1:length(Ind2))
    #   if(Ind2[i]!=ind[i]) {
    #     ErrorMSG<<-paste("Input data series not continuous at:",Ind2[i],ind[i],"\n")
    #     return(-1)
    #   }
    IY0<-ind[is.na(itable[,4])==F]
    IY0flg<-rep(0,length(IY0))
    Y0<-itable[is.na(itable[,4])==F,4]
    olflg<-is.na(itable[,4])==F
    Iyr<-floor(IY0/10000)
    Imd<-IY0-Iyr*10000
    Ti<-IY0
    for(i in 1:length(IY0)){
        ith<-match(Imd[i],Icy)
        Ti[i]<-(Iyr[i]-iyrbegin)*Nt+ith
    }
    for(i in 1:(length(IY0)-1)){
        if(Ti[i+1]-Ti[i]==1) IY0flg[i]<-1
    }
    if(sum(IY0flg)<10){  # too few available data for autocorlh
        ErrorMSG<<-paste("Too many missing values in ", "to estimate autocorrelation\n")
        return(-1)
    }
    itable <- itable[is.na(itable[, 4]) == F, ]
    assign("ori.itable", ori.itable, envir = .GlobalEnv)
    assign("itable", itable, envir = .GlobalEnv)
    assign("Ti", Ti, envir = .GlobalEnv) # Time index for LS fitting
    assign("Y0", Y0, envir = .GlobalEnv) # Data series for Base
    assign("IY0", IY0, envir = .GlobalEnv) # Cycle index for Base
    assign("Imd", Imd, envir = .GlobalEnv) # Cycle index for Base
    assign("IY0flg", IY0flg, envir = .GlobalEnv) # continuous flag for Base
    assign("Icy", Icy, envir = .GlobalEnv) # Cycle index
    assign("Nt", Nt, envir = .GlobalEnv) # Cycle length
    assign("olflg", olflg, envir = .GlobalEnv)
    return(0)
}
