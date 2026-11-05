#' read input of RHtestsV4
#'
#' Remove NA values and flag `y[t]` as 0 when it has not corresponding `y[t+1]` values.
#'
#' @inheritParams FindU
#'
#' @return
#' An S7 `RHState` object containing the series, seasonal cycle and critical values.
#' `Nt` is 1, 12 or 365 for annual, monthly or daily series, respectively.
#' @examples
#' infile <- system.file("extdata/Example1.dat", package = "RHtests")
#' Read(infile, "-999.99")
#' @export
Read <- function(InSeries, MissingValueCode = "-999.99", plev = 0.95) {
  if (is.null(InSeries)) stop("InSeries is required", call. = FALSE)

  itable <- if (is.character(InSeries)) {
    sep <- ifelse(is.csv(InSeries), ",", " ")
    fread(InSeries, sep = sep, na.strings = MissingValueCode)
  } else {
    InSeries
  }
  itable %<>% set_colnames(c("year", "month", "day", "data")) %>% as.data.table()

  # keep input base data as ori.itable
  ori.itable <- itable
  isnot_0229 <- itable[, !(month == 2 & day == 29)]
  ooflg <- itable[, which(isnot_0229 & !is.na(data))]

  # get rid of Feb 29th data

  itable <- itable[isnot_0229, ]
  iyrbegin <- itable$year[1]
  imdbegin <- itable$month[1] * 100 + itable$day[1]
  n <- nrow(itable)
  iyrend <- itable$year[n]
  imdend <- itable$month[n] * 100 + itable$day[n]

  # check input data (both base and ref), no jump with begin and end
  Icy <- itable[, sort(unique(month * 100 + day))]
  # Icy  <- sort(unique(ori.itable[ooflg, 2]*100+ori.itable[ooflg,3]))
  Ind2 <- iyrbegin * 10000 + Icy[Icy >= imdbegin] # first year

  if (iyrend > (iyrbegin + 1)) {
    for (i in (iyrbegin + 1):(iyrend - 1)) {
      Ind2 <- c(Ind2, i * 10000 + Icy)
    }
  }
  Ind2 <- c(Ind2, iyrend * 10000 + Icy[Icy <= imdend])
  Nt <- length(Icy)

  ind <- ori.itable[ooflg, year * 10000 + month * 100 + day]
  if (sum(!(ind %in% Ind2)) > 0 || all.equal(ind, sort(ind)) != TRUE) {
    stop("Input data series not continuous", call. = FALSE)
  }
  # 允许缺测，但有效观测日期须按时间顺序排列。
  IY0 <- ind
  IY0flg <- rep(0, length(IY0))
  Y0 <- itable[!is.na(data), data]
  Iyr <- floor(IY0 / 10000)
  Imd <- IY0 - Iyr * 10000
  Ti <- IY0
  for (i in 1:length(IY0)) {
    ith <- match(Imd[i], Icy)
    Ti[i] <- (Iyr[i] - iyrbegin) * Nt + ith
  }

  for (i in 1:(length(IY0) - 1)) {
    if (Ti[i + 1] - Ti[i] == 1) IY0flg[i] <- 1
  }
  if (sum(IY0flg) < 10) { # too few available data for autocorlh
    stop("Too many missing values to estimate autocorrelation", call. = FALSE)
  }

  N <- length(Y0)
  critical_values <- readPFtable(N, plev)

  itable <- itable[!is.na(data), ]
  state <- listk(
    ori.itable = as.data.frame(ori.itable), itable = as.data.frame(itable),
    ooflg, Ti, Y0, IY0, Imd, IY0flg, Icy, Nt
  )
  new_rh_state(state, critical_values, plev = plev)
}

is.csv <- function(names) {
  nlen <- nchar(names)
  if (substr(names, nlen - 2, nlen) == "csv" || substr(names, nlen - 2, nlen) == "CSV") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


ReadDLY.g <- function(InSeries, MissingValueCode) {
  if (!file.exists(InSeries)) stop("Input datafile does not exist: ", InSeries, call. = FALSE)
  if (is.csv(InSeries)) {
    itmp <- try(read.table(InSeries,
      sep = ",", header = FALSE, na.strings = MissingValueCode,
      colClasses = rep("real", 4)
    ), silent = TRUE)
    if (inherits(itmp, "try-error")) {
      stop(as.character(itmp), call. = FALSE)
    } else {
      itable <- itmp
    }
  } else {
    itmp <- try(read.table(InSeries,
      sep = "", header = FALSE, na.strings = MissingValueCode,
      colClasses = rep("real", 4)
    ), silent = TRUE)
    if (inherits(itmp, "try-error")) {
      stop(as.character(itmp), call. = FALSE)
    } else {
      itable <- itmp
    }
  }
  if (ncol(itable) != 4) stop(InSeries, " must have 4 columns", call. = FALSE)
  colnames(itable) <- c("id1", "id2", "id3", "data")

  iyrbegin <- itable[1, 1]
  imdbegin <- itable[1, 2] * 100 + itable[1, 3]
  iyrend <- itable[dim(itable)[1], 1]
  imdend <- itable[dim(itable)[1], 2] * 100 + itable[dim(itable)[1], 3]
  # keep input base data as ori.itable
  ori.itable <- itable
  # check input data (both base and ref), no jump with begin and end
  Icy <- sort(unique(itable[, 2] * 100 + itable[, 3]))
  Ind2 <- iyrbegin * 10000 + Icy[Icy >= imdbegin] # first year
  # if(iyrend>(iyrbegin+1)) for(i in (iyrbegin+1):(iyrend-1))
  #   Ind2<-c(Ind2,i*10000+Icy)
  # Ind2<-c(Ind2,iyrend*10000+Icy[Icy<=imdend])
  Nt <- length(Icy)
  Nall <- dim(itable)[1]
  ind <- itable[, 1] * 10000 + itable[, 2] * 100 + itable[, 3]
  # 按有效观测日期构造拟合索引。
  IY0 <- ind[is.na(itable[, 4]) == FALSE]
  IY0flg <- rep(0, length(IY0))
  Y0 <- itable[is.na(itable[, 4]) == FALSE, 4]
  olflg <- is.na(itable[, 4]) == FALSE
  Iyr <- floor(IY0 / 10000)
  Imd <- IY0 - Iyr * 10000
  Ti <- IY0
  for (i in 1:length(IY0)) {
    ith <- match(Imd[i], Icy)
    Ti[i] <- (Iyr[i] - iyrbegin) * Nt + ith
  }
  for (i in 1:(length(IY0) - 1)) {
    if (Ti[i + 1] - Ti[i] == 1) IY0flg[i] <- 1
  }
  if (sum(IY0flg) < 10) { # too few available data for autocorlh
    stop("Too many missing values to estimate autocorrelation", call. = FALSE)
  }
  itable <- itable[is.na(itable[, 4]) == FALSE, ]
  listk(ori.itable, itable, Ti, Y0, IY0, Imd, IY0flg, Icy, Nt, olflg)
}
