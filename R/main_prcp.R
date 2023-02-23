BCtrans <- function(P, lambda) {
  if (sum(P <= 0) > 0) stop("non-positive value encountered, y<=0!")
  y <- if (lambda == 0) log(P) else (P**lambda - 1) / lambda
  return(y)
}

IVBCtrans <- function(P, lambda) {
  y <- if (lambda == 0) exp(P) else (P * lambda + 1)**(1 / lambda)
  return(y)
}

getMtrendFdly <- function(idata) {
  cmon <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Agu", "Sep", "Oct", "Nov", "Dec")
  if (dim(idata)[2] != 4) stop("input data does not contain 4 columns in getMtrendFdly())")
  colnames(idata) <- c("year", "month", "day", "data")
  yrs <- unique(idata[, "year"])
  mmean <- rep(NA, 12)
  tdata <- NULL
  for (yth in 1:length(yrs)) {
    itmp <- ori.itable[ori.itable[, 1] == yrs[yth], ]
    for (mon in 1:12) {
      mp <- NA
      mtmp <- itmp[itmp[, 2] == mon, ]
      if (length(mtmp) > 0) {
        if (dim(mtmp)[1] > 20 & sum(is.na(mtmp[, 4])) <= 3) mp <- sum(idata[idata[, 1] == yrs[yth] & idata[, 2] == mon, 4])
      }
      tdata <- rbind(tdata, c(yrs[yth], mon, mp))
    }
  }
  t1data <- tdata
  for (i in 1:12) {
    if (sum(is.na(tdata[tdata[, 2] == i, 3]) == F) < 5) print(paste("monthly data too few at", cmon[i]))
    mmean[i] <- mean(tdata[tdata[, 2] == i, 3], na.rm = T)
    t1data[t1data[, 2] == i, 3] <- t1data[t1data[, 2] == i, 3] - mmean[i]
  }
  y <- t1data[, 3]
  ind <- which(!is.na(y))
  y <- y[ind]
  # print(y)
  summary(lm(y ~ ind))$coef[2, 1]
}
