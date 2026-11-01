# 搜索参考序列检验中的候选突变点。
search_pmt_candidates <- function(Y0, Ips, Ns, Nmin, phi, PTmax, plev) {
  # take cor==0, search for all possible changepoints
  tt <- TRUE
  while (tt) {
    tt <- FALSE
    Ips0 <- NULL
    for (Iseg in 1:(Ns + 1)) {
      I0 <- if (Iseg == 1) 0 else Ips[Iseg - 1]
      I2 <- Ips[Iseg]
      otmp <- PTKI0I2(Y0, I0, I2, Nmin)
      if (otmp$prob > 0) Ips0 <- sort(c(Ips0, otmp$Ic))
    }
    # finished find new possible changepoints, sorted in Ips0
    tt1 <- TRUE
    while (tt1) { # check Ips0 and insert new break points
      if (length(Ips0) == 0) {
        tt1 <- FALSE
      } else {
        PTx.mx <- (-9999)
        PTx95.mx <- .001
        prob.mx <- 0
        for (i in 1:length(Ips0)) { # search Ips0 series, find max prob
          Ips1 <- sort(c(Ips, Ips0[i]))
          Ic <- Ips0[i]
          id <- match(Ic, Ips1)
          if (id == length(Ips1)) {
            print(Ips1)
            print(id)
            stop("error in FindSteps")
          }
          I0 <- if (id == 1) 1 else Ips1[id - 1] + 1
          I2 <- Ips1[id + 1]
          Nseg <- I2 - I0
          PTx95 <- getPTx95(0, (Nseg - 1), phi, PTmax)
          Pk1 <- Pk.PMT(Nseg)
          otmp <- PTKIc(Y0[I0:I2], Pk1, Ic - I0 + 1)
          if (otmp$prob < plev) {
            Ips0[i] <- 0
          } else if (otmp$PTk / PTx95 > PTx.mx / PTx95.mx) {
            PTx.mx <- otmp$PTk
            PTx95.mx <- PTx95
            Imx <- Ic
            inc <- i
            prob.mx <- otmp$prob
          }
        }
        if (prob.mx >= plev) {
          Ips <- sort(c(Ips, Imx)) # insert new point into Ips
          Ns <- Ns + 1
          Ips0 <- Ips0[-inc] # exclude co-responding point in Ips0
          tt <- TRUE
        } else {
          tt1 <- FALSE
        } # finish inserting new points into Ips
        Ips0 <- Ips0[Ips0 != 0]
      }
    }
  } # end of searching for all possible changepoints
  list(Ips = Ips, Ns = Ns)
}
