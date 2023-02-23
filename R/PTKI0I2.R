PTKI0I2 <- function(Y0, I0, I2, Nmin) {
  # search new breakpoint of Y0[(I0+1):I2] using function PTK()
  # output: Ic -- breakpoint, prob and PTx
  Y <- Y0[(I0 + 1):I2]
  N <- length(Y)
  oout <- list(prob = -1, Ic = I0, PTx = -9999.9)

  if (N >= (Nmin * 2)) {
    Pk0 <- Pk.PMT(N)
    otmp <- PTK(Y, Pk0, Nmin)
    oout$Ic <- I0 + otmp$KPx
    oout$prob <- pt(otmp$Tx, (N - 2))
    oout$PTx <- otmp$PTx
  }
  return(oout)
}

PTK <- function(Y, Pk, Nmin) {
  #  search input vector, return PTx.max and corresponding KPx
  PTx <- (-9999.9)
  KPx <- 0
  N <- length(Y)
  ALL <- sum(Y)

  for (k in Nmin:(N - Nmin)) {
    EY1 <- mean(Y[1:k])
    EY2 <- (ALL - EY1*k) / (N - k) # mean(Y[(k + 1):N])
    var <- sum(c((Y[1:k] - EY1)^2, (Y[(k + 1):N] - EY2)^2))
    std <- sqrt(var / (N - 2))
    Tk <- sqrt(k * (N - k) / N) * abs(EY1 - EY2) / std
    PTk <- Tk * Pk[k]
  
    # Ipaper::fprintf("var = %f, sd = %f, Tk = %f, PTk = %f \n", 
    #   var, std, Tk, PTk)
    if (PTk > PTx) {
      PTx <- PTk
      KPx <- k
      Tx <- Tk
    }
  }
  list(PTx = PTx, KPx = KPx, Tx = Tx)
}
