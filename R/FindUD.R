#' FindUD
#' 
#' @inheritParams FindU
#' 
#' @param InCs A data.frame of detected turning points returned from [FindU],
#' with the columns of `c("kind", "Idc", "date", "probL", "probU", "plev",
#' "PFx", "PFx95l", "PFx95h")`.
#' 
#' @examples 
#' f = system.file("extdata/Example1.dat", package = "RHtests")
#' r1 = FindU(f)
#' 
#' r2 = FindUD(InCs = r1$TP, output = "OUTPUT/example1")
#' print(r2)
#' 
#' plot_output(r2$fit)
#' @export
FindUD <- function(
    InSeries = NULL, InCs, output, MissingValueCode = "-999.99",
    GUI = FALSE, plev = 0.95, Iadj = 10000, Mq = 10, Ny4a = 0, is_plot = FALSE) {
  
  if (!is.null(InSeries)) data <- Read(InSeries, MissingValueCode) # data not used
  # all necessary variables have been in global env
  
  Debug <- TRUE
  ErrorMSG <- NA
  assign("ErrorMSG", ErrorMSG, envir = .GlobalEnv)
  flog <- paste(output, ".log", sep = "")
  Nmin <- 10
  if (Ny4a > 0 & Ny4a <= 5) Ny4a <- 5

  assign("Nmin", Nmin, envir = .GlobalEnv)
  N <- length(Y0)
  Nadj <- Ny4a * Nt

  # change InCs
  TP <- if (is.character(InCs)) fread(InCs) else InCs

  if (is.Date(TP$date)) TP$date %<>% date2num()
  Ns <- nrow(TP) # number of changing points
  if (is.null(Ns)) Ns <- 0

  Pk0 <- Pk.PMFT(N)
  oout <- rmCycle(itable)
  Y1 <- oout$Base
  EB <- oout$EB
  assign("EB", EB, envir = .GlobalEnv)
  if (length(EB) != length(Icy)) {
    ErrorMSG <<- paste(
      "Annual cycle length (from non-missing data) differ from original dataset",
      "\n", get("ErrorMSG", env = .GlobalEnv), "\n")
    if (!GUI) cat(ErrorMSG)
    return(-1)
  }

  ofileIout <- paste(output, "_pCs.txt", sep = "")
  ofileMout <- paste(output, "_mCs.txt", sep = "")

  ofilePdf <- paste(output, "_UD.pdf", sep = "")
  ofileSout <- paste(output, "_UDstat.txt", sep = "")
  file.create(ofileIout)
  file.create(ofilePdf)
  file.create(ofileSout)
  cat(paste("The nominal level of confidence (1-alpha)=", plev, "\n"), file = ofileSout)
  cat(paste("Input data filename:", "N=", N, "\n"), file = ofileSout, append = T)

  Ip0 <- N
  oout <- LSmultiRedCycle(Y1, Ti, Ip0, 1)
  beta0 <- oout$trend
  betaL0 <- oout$betaL
  betaU0 <- oout$betaU
  meanhat0 <- oout$meanhat
  Ehat0 <- mean(meanhat0)
  p.tr0 <- oout$p.tr
  corD <- oout$cor
  corDL <- oout$corl
  corDU <- oout$corh

  if (Ns == 0) { # no input changepoints no input changepoints
    oout <- PMFT(Y1, Ti, Pk0)
    I0 <- 0
    I2 <- oout$KPx
    I4 <- N

    oout1 <- PMFxKxI0I2(Y1, Ti, I0, I2)
    I1 <- oout1$Ic
    oout2 <- PMFxKxI0I2(Y1, Ti, I2, I4)
    I3 <- oout2$Ic
    oout3 <- PMFxKxI0I2(Y1, Ti, I1, I3)
    I2 <- oout3$Ic

    Ns <- 1
    Ips <- c(I1, N)
    if (I1 > 0) {
      otmp <- LSmultiple(Y1, Ti, Ips)
      resi <- otmp$resi
      fitted <- otmp$fitted
      otmp <- Rphi(resi, Ips, Ns)
      cor1 <- otmp$cor
      corL1 <- otmp$corl
      corU1 <- otmp$corh
      W <- otmp$W + fitted
      otmp <- PMFxKc(Y1, Ti, I0, I4, I1)
      PFx1 <- otmp$PFc
      otmp <- PMFxKc(W, Ti, I0, I4, I1)
      prob1 <- otmp$prob
    } else {
      prob1 <- 0
      PFx1 <- 0
    }

    Ips <- c(I2, N)
    if (I2 > 0) {
      otmp <- LSmultiple(Y1, Ti, Ips)
      resi <- otmp$resi
      fitted <- otmp$fitted
      otmp <- Rphi(resi, Ips, Ns)
      cor2 <- otmp$cor
      corL2 <- otmp$corl
      corU2 <- otmp$corh
      W <- otmp$W + fitted
      otmp <- PMFxKc(Y1, Ti, I0, I4, I2)
      PFx2 <- otmp$PFc
      otmp <- PMFxKc(W, Ti, I0, I4, I2)
      prob2 <- otmp$prob
    } else {
      prob2 <- 0
      PFx2 <- 0
    }

    Ips <- c(I3, N)
    if (I3 > 0) {
      otmp <- LSmultiple(Y1, Ti, Ips)
      resi <- otmp$resi
      fitted <- otmp$fitted
      otmp <- Rphi(resi, Ips, Ns)
      cor3 <- otmp$cor
      corL3 <- otmp$corl
      corU3 <- otmp$corh
      W <- otmp$W + fitted
      otmp <- PMFxKc(Y1, Ti, I0, I4, I3)
      PFx3 <- otmp$PFc
      otmp <- PMFxKc(W, Ti, I0, I4, I3)
      prob3 <- otmp$prob
    } else {
      prob3 <- 0
      PFx3 <- 0
    }

    tmp <- sort(c(PFx1, PFx2, PFx3), decreasing = T, index.return = T)
    PFx.mx <- tmp$x[1]
    prob.mx <- c(prob1, prob2, prob3)[tmp$ix[1]]
    Imx <- c(I1, I2, I3)[tmp$ix[1]]
    if (prob.mx < plev) {
      #     cat("PMF finds the series to be homogeneous!\n",file=ofileIout)
      cat(paste(0, "changepoints in Series", "\n"), file = ofileIout)
      cat("PMF finds the series to be homogeneous!\n")
      #     return()
      Ns <- 0
      Ips <- N
      Ids <- c(0)
    } else {
      Ns <- 1
      Ips <- c(Imx, N)
      Ids <- c(0, 1)
    }
  } else {
    Ips <- c(rep(0, Ns), N) #
    Ids <- rep(0, Ns)

    for (i in 1:Ns) { # using YYYYMMDD as index, searching for the largest
      # date less or equal to given YYYYMMDD
      # ymdtmp <- as.numeric(substr(itmp[i+1],7,16))
      ymdtmp <- TP$date[i] # %>% format("%Y%m%d") %>% as.numeric()
      it <- match(ymdtmp, IY0)
      if (!is.na(it)) {
        Ips[i] <- it
      } else {
        Ips[i] <- max(c(1:N)[IY0 <= ymdtmp])
      }
      # Ids[i] <- as.numeric(substr(itmp[i+1],1,1))
    }
    Ids <- TP$kind

    if (sum(is.na(Ips)) > 0 | !identical(Ips, sort(Ips))) {
      ErrorMSG <<- paste("FindUD: Ips read in from ", InCs, "error:")
      for (i in 1:Ns) {
        ErrorMSG <<- paste(get("ErrorMSG", env = .GlobalEnv), Ips[i])
      }
      ErrorMSG <<- paste(get("ErrorMSG", env = .GlobalEnv), "\n")
      if (!GUI) cat(ErrorMSG)
      return(-1)
    }
  }

  if (Ns > 0) {
    Ips.i <- Ips
    Niter <- 0
    # start search for all possible changepoints
    tt <- TRUE
    while (tt) {
      Niter <- Niter + 1
      tt <- FALSE
      Ips0 <- NULL
      for (i in 1:(Ns + 1)) {
        I0 <- if (i == 1) 0 else Ips[i - 1]
        I2 <- Ips[i]
        otmp <- PMFxKxI0I2(Y1, Ti, I0, I2)
        if (otmp$prob > 0) Ips0 <- sort(c(Ips0, otmp$Ic))
      }
      # estimate p-value of each changepoint in series Ips0, find the most significant
      # changepoint Iseg.mx
      tt1 <- TRUE
      while (tt1) {
        if (length(Ips0) == 0) {
          tt1 <- FALSE
        } else {
          Iseg.mx <- 0
          prob.mx <- (-1)
          probL.mx <- (-1)
          PFx.mx <- (-1)
          for (i in 1:length(Ips0)) {
            Ips1 <- sort(c(Ips, Ips0[i]))
            ith <- match(Ips0[i], Ips1)
            otmp <- PMFxIseg(Y1, Ti, Ips1, ith)
            probL <- min(c(otmp$probL, otmp$probU, otmp$prob))
            probU <- max(c(otmp$probL, otmp$probU, otmp$prob))
            PFx <- otmp$PFx
            if (probU < plev) {
              Ips0[i] <- 0
            } else if (PFx > PFx.mx) {
              prob.mx <- otmp$prob
              probL.mx <- probL
              Iseg.mx <- Ips0[i]
              PFx.mx <- PFx
            }
          }
          if (probL.mx >= plev) {
            Ips <- sort(c(Ips, Iseg.mx))
            Ns <- Ns + 1
            Ips0 <- Ips0[Ips0 != Iseg.mx]
            tt <- TRUE
          } else {
            tt1 <- FALSE
          }
          Ips0 <- Ips0[Ips0 != 0]
        }
      }
    }
    Ids0 <- rep(NA, length(Ips))
    for (i in 1:length(Ips)) {
      if (Ips[i] %in% Ips.i) {
        Ids0[i] <- Ids[Ips.i == Ips[i]]
      } else {
        Ids0[i] <- 0
      }
    }
    Ids <- Ids0
    # Ids<-as.integer(Ips%in%Ips.i)
    tt <- TRUE
    while (tt) {
      tt <- FALSE
      probL.mn <- 9999
      Iseg.mn <- 0
      for (i in 1:Ns) {
        if (Ids[i] == 0) { # check those un-documented
          Ips0 <- Ips[-i]
          otmp <- PMFxIseg(Y1, Ti, Ips, i)
          probL <- min(otmp$probL, otmp$probU)
          if (probL < probL.mn) {
            Iseg.mn <- i
            probL.mn <- probL
          }
        } # end if documented
      }
      if (Iseg.mn > 0 & probL.mn < plev) {
        Ips <- Ips[-Iseg.mn]
        Ids <- Ids[-Iseg.mn]
        Ns <- Ns - 1
        if (Ns > 0) tt <- TRUE
      }
    }
  }
  # all changepoints are significant, calculate stats and output

  if (Ns > 0) {
    Nsegs <- Ips - c(0, Ips[1:Ns])
    Iseg.longest <- sort(Nsegs, index = T, decreasing = T)$ix[1]
  } else {
    Iseg.longest <- 0
  }

  if (Iadj > (Ns + 1) | Iseg.longest == 0) {
    Iseg.adj <- Ns + 1
  } else if (Iadj == 0) {
    Iseg.adj <- Iseg.longest
  } else {
    Iseg.adj <- Iadj
  }

  otmp <- LSmultiRedCycle(Y1, Ti, Ips, Iseg.adj)
  Y1 <- otmp$Y0
  cor <- otmp$cor
  corl <- otmp$corl
  corh <- otmp$corh
  df <- (N - 2 - Nt - Ns)
  pcor <- pt(abs(cor) * sqrt(df / (1 - cor^2)), df)
  Rf <- otmp$resi
  W <- otmp$W
  WL <- otmp$WL
  WU <- otmp$WU
  EB1 <- otmp$EB

  itmp1 <- cbind(EB1, Icy)
  itmp2 <- cbind(1:N, Imd)
  colnames(itmp2) <- c("idx", "Icy")

  itmp <- merge(itmp1, itmp2, by = "Icy")
  EBfull <- itmp[order(itmp[, "idx"]), "EB1"]
  EEB <- mean(EB1, na.rm = T)

  if (Ns > 0) {
    Rb <- Y1 - otmp$trend * Ti + EBfull
    QMout <- QMadjGaussian(Rb, Ips, Mq, Iseg.adj, Nadj)
    B <- QMout$PA
    cat(paste("Nseg_shortest =", QMout$Nseg.mn, "; Mq = ", QMout$Mq, "; Ny4a = ", Ny4a, "\n"),
      file = ofileSout, append = T
    )
    cat(paste(
      "\n Adjust to segment", Iseg.adj, ": from",
      if (Iseg.adj == 1) 1 else Ips[Iseg.adj - 1] + 1,
      "to", Ips[Iseg.adj], "\n"
    ), file = ofileSout, append = T)
    #   cat("#Fcat, DP (CDF and Differnces in category mean)\n",file=ofileSout,
    #       append=T)
    if (QMout$Mq > 1) {
      oline <- paste("#Fcat: frequency category boundaries\n",
        "#DP: Difference in the category means\n#", sep = "")
      for (i in 1:Ns) oline <- paste(oline, paste("Fcat.Seg", i, sep = ""), 
        paste("DP.Seg", i, sep = ""))
      
      oline <- paste(oline, "\n")
      cat(oline, file = ofileSout, append = T)

      write.table(round(QMout$osmean, 4),
        file = ofileSout, append = T, row.names = F, col.names = F)
      
      for (i in 1:(Ns + 1)) {
        I1 <- if (i == 1) 1 else Ips[i - 1] + 1
        I2 <- Ips[i]
        if (i != Iseg.adj) {
          cat(paste("Seg. ", i, ": mean of QM-adjustments =", round(mean(QMout$W[I1:I2]), 4),
            "\n", sep = ""), file = ofileSout, append = T)
        }
      }
    }
  } else {
    B <- Y1 - otmp$trend * Ti + EBfull
  }

  adj <- otmp$Y0 + EBfull
  # B<-B+EBfull+otmp$trend*Ti
  B <- B + otmp$trend * Ti

  cat(file = ofileSout, paste(
    " Ignore changepoints -> trend0 =",
    round(beta0, 6), "(", round(betaL0, 6), ",", round(betaU0, 6), ") (p=",
    round(p.tr0, 4), "); cor=", round(corD, 4), "(", round(corDL, 4),
    ",", round(corDU, 4), ")\n"
  ), append = T)
  cat("Common trend TPR fit to the seasonal-cycle-adjusted Base series:\n",
    file = ofileSout, append = T
  )
  cat(
    paste(
      "#steps= ", Ns, "; trend=", round(otmp$trend, 6), "(",
      round(otmp$betaL, 6), ",", round(otmp$betaU, 6), ") (p=",
      round(otmp$p.tr, 4), "); cor=",
      round(cor, 4), "(", round(corl, 4), ",", round(corh, 4), ")",
      round(pcor, 4), "\n"
    ),
    file = ofileSout, append = T
  )
  if (Ns > 0) {
    for (i in 1:(Ns + 1)) {
      I1 <- if (i == 1) 1 else Ips[i - 1] + 1
      I2 <- Ips[i]
      Delta <- otmp$mu[Iseg.adj] - otmp$mu[i]
      adj[I1:I2] <- adj[I1:I2] + Delta
      stepsize <- otmp$mu[i + 1] - otmp$mu[i]
      cat(paste(Ips[i], IY0[Ips[i]], "stepsize=", round(stepsize, 4), "\n"),
        file = ofileSout, append = T
      )
    }
  }

  oR <- Y1 - otmp$meanhat
  oR[2:N] <- oR[2:N] - oR[1:(N - 1)] * cor
  Ehat <- mean(otmp$meanhat)
  meanhat0 <- meanhat0 - Ehat0 + Ehat

  if (is_plot) {
    plot_FindUD(
      oout, ofilePdf, EBfull, EEB, B, QMout, Ms, Mq, Ns, adj,
      Ips, Iseg.adj, otmp
    )
  }

  odata <- matrix(NA, dim(ori.itable)[1], 10) %>% set_colnames(fitdata_varnames_noref)
  # odata[ooflg,1]<-Ti
  odata[, 1] <- c(1:dim(ori.itable)[1])
  odata[, 2] <- ori.itable[, 1] * 10000 + ori.itable[, 2] * 100 + ori.itable[, 3]
  # odata[ooflg,3] <- round(otmp$Y0+EBfull,4)
  odata[, 3] <- ori.itable[, 4]
  odata[ooflg, 4] <- round(otmp$meanhat + EEB, 4)
  odata[ooflg, 5] <- round(adj, 4)
  odata[ooflg, 6] <- round(otmp$Y0, 4)
  odata[ooflg, 7] <- round(otmp$meanhat, 4)
  odata[ooflg, 8] <- round(otmp$meanhat + EBfull, 4)
  if (Ns > 0) if (QMout$Mq > 1) odata[ooflg, 9] <- round(B, 4)
  # odata[ooflg,9]<-round(oR,4)
  # odata[ooflg,10]<-round(Rb,4)
  odata[ooflg, 10] <- round(meanhat0, 4)

  Imd1 <- ori.itable[, 2] * 100 + ori.itable[, 3]
  if (sum(is.na(ori.itable[, 4]) == F & Imd1 == 229) > 0) {
    if (Ns > 0) {
      tdata <- ori.itable[is.na(ori.itable[, 4]) == F, ]
      IY1 <- tdata[, 1] * 10000 + tdata[, 2] * 100 + tdata[, 3]
      Ips.ymd <- IY0[Ips]
      Ips.1 <- rep(NA, Ns + 1)
      for (i in 1:Ns) Ips.1[i] <- c(1:length(IY1))[IY1 == Ips.ymd[i]]
      Ips.1[Ns + 1] <- length(IY1)
      #     Ips.1<-c(1:length(IY1))[Ips.ymd==IY1]
      Imd2 <- tdata[, 2] * 100 + tdata[, 3]
      Ids.leap <- c(1:length(Imd2))[Imd2 == 229]
      Nl <- length(Ids.leap)
      Rb <- Y1 - otmp$trend * Ti + EBfull
      Rb1 <- tdata[, 4]
      Rb1[-Ids.leap] <- Rb
      Ti1 <- rep(NA, length(IY1))
      Ti1[-Ids.leap] <- Ti
      for (i in 1:length(Ids.leap)) {
        Rb1[Ids.leap[i]] <- tdata[Ids.leap[i], 4] + Rb1[Ids.leap[i] - 1] - tdata[Ids.leap[i] - 1, 4]
        Ti1[Ids.leap[i]] <- Ti1[Ids.leap[i] - 1]
      }
      if (QMout$Mq > 1) {
        B1 <- QMadjGaussian(Rb1, Ips.1, Mq, Iseg.adj, Nadj)$PA
        B1 <- B1 + otmp$trend * Ti1
        B1.leap <- B1[Ids.leap]
        odata[is.na(odata[, 3]) == F & Imd1 == 229, 9] <- round(B1.leap, 4)
      }
    } else {
      odata[Imd1 == 229, 9] <- odata[Imd1 == 229, 3]
    }
    Ids.leapo <- c(1:dim(ori.itable)[1])[is.na(ori.itable[, 4]) == F & Imd1 == 229]
    for (jth in 1:length(Ids.leapo)) {
      kth <- Ids.leapo[jth]
      if (Ns > 0) {
        k1th <- if (odata[kth - 1, 2] %in% IY0[Ips]) (kth + 1) else (kth - 1)
      } else {
        k1th <- kth - 1
      }
      for (pth in c(4, 7, 8, 10)) odata[kth, pth] <- odata[k1th, pth]
      for (pth in c(5, 6)) {
        delta1 <- odata[k1th, 3] - odata[k1th, pth]
        odata[kth, pth] <- odata[kth, 3] - delta1
      }
    }
  }

  RW <- LSmultiple(W, Ti, Ips)$resi
  RWL <- LSmultiple(WL, Ti, Ips)$resi
  RWU <- LSmultiple(WU, Ti, Ips)$resi

  d_TP <- list()
  if (Ns > 0) {
    # d_TP = foreach(i = 1:Ns) %do% {
    for (i in 1:Ns) {
      I1 <- if (i == 1) 1 else Ips[i - 1] + 1
      I3 <- Ips[i + 1]
      Ic <- Ips[i]
      Id <- Ids[i]
      Nseg <- I3 - I1 + 1

      PFx95 <- getPFx95(cor, Nseg)
      PFx95l <- getPFx95(corl, Nseg)
      PFx95h <- getPFx95(corh, Nseg)
      SSEf.Iseg <- sum(Rf[I1:I3]^2)
      Ips0 <- Ips[-i]

      otmp <- LSmultiple(Y1, Ti, Ips0)
      SSE0.Iseg <- sum(otmp$resi[I1:I3]^2)
      Fx <- (SSE0.Iseg - SSEf.Iseg) * (Nseg - 3) / SSEf.Iseg

      Pk0 <- Pk.PMFT(Nseg)
      PFx <- Fx * Pk0[Ic - I1 + 1]

      otmp <- LSmultiple(W, Ti, Ips0)
      SSE0.Iseg <- sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg <- sum(RW[I1:I3]^2)
      Fx <- (SSE0.Iseg - SSEf.Iseg) * (Nseg - 3) / SSEf.Iseg

      if (Fx <= 0) {
        PFx <- 0
        Fx <- 0
        prob <- 0
      } else {
        prob <- pf(Fx, 1, Nseg - 3)
      }

      otmp <- LSmultiple(WL, Ti, Ips0)
      SSE0.Iseg <- sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg <- sum(RWL[I1:I3]^2)
      Fx <- (SSE0.Iseg - SSEf.Iseg) * (Nseg - 3) / SSEf.Iseg
      probL0 <- if (Fx < 0) 0 else pf(Fx, 1, Nseg - 3)

      otmp <- LSmultiple(WU, Ti, Ips0)
      SSE0.Iseg <- sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg <- sum(RWU[I1:I3]^2)
      Fx <- (SSE0.Iseg - SSEf.Iseg) * (Nseg - 3) / SSEf.Iseg
      probU0 <- if (Fx < 0) 0 else pf(Fx, 1, Nseg - 3)

      probL <- min(probL0, probU0)
      probU <- max(probL0, probU0)

      ## type_TP的判断标准
      # browser()
      Idc <- guess_sign_level(Id, plev, probL, probU, PFx, PFx95l, PFx95h)

      cat(
        paste(sprintf("%1.0f", as.numeric(Id)), " ",
          sprintf("%-4.4s", Idc),
          sprintf("%10.0f", IY0[Ic]), " (",
          sprintf("%6.4f", probL), "-",
          sprintf("%6.4f", probU), ")",
          sprintf("%6.3f", plev),
          sprintf("%10.4f", PFx), " (",
          sprintf("%10.4f", PFx95l), "-",
          sprintf("%10.4f", PFx95h), ")\n",
          sep = ""
        ),
        file = ofileIout,
        append = TRUE
      )
      cat(
        paste("PMF : c=", sprintf("%4.0f", Ic),
          "; (Time ", sprintf("%10.0f", IY0[Ic]),
          "); Type= ", sprintf("%4.0f", as.numeric(Id)),
          "; p=", sprintf("%10.4f", prob),
          "(", sprintf("%10.4f", probL),
          "-", sprintf("%10.4f", probU),
          "); PFmax=", sprintf("%10.4f", PFx),
          "; CV95=", sprintf("%10.4f", PFx95),
          "(", sprintf("%10.4f", PFx95l),
          "-", sprintf("%10.4f", PFx95h),
          "); Nseg=", sprintf("%4.0f", Nseg), "\n",
          sep = ""
        ),
        file = ofileSout, append = T
      )
      d_TP[[i]] <- data.table(kind = Id, Idc, date = IY0[Ic], probL, probU, plev, PFx, PFx95l, PFx95h)
    }

    if (Ns == 0) cat("PMF finds the series to be homogeneous!\n")
    cat(paste("# ", Ns, "changepoints in Series", "\n"), file = ofileIout)

    if (!is_empty(d_TP)) {
      d_TP %<>% do.call(rbind, .)
      d_TP[, 4:9] %<>% lapply(round, digits = 4)
      fwrite(d_TP, ofileIout, append = TRUE, col.names = TRUE)
    }
  }

  odata %<>% as.data.table()
  odata$date %<>% num2date()
  write.table(file = paste0(output, "_UD.dat"), odata, na = MissingValueCode, col.names = TRUE, row.names = F)

  if (GUI) {
    return(0)
  } else {
    file.copy(from = ofileIout, to = ofileMout, overwrite = TRUE)
    # cat("FindUD finished successfully...\n")
    list(fit = odata, TP = d_TP) # fit = odata,
  }
}

guess_sign_level <- function(Id, plev, probL, probU, PFx, PFx95l, PFx95h) {
  if (Id == 0) { # type-0 changepoints
    if (probU < plev) {
      Idc <- "No  "
    } else if (probL < plev & probU >= plev) {
      Idc <- "?   "
    } else if (probL >= plev) Idc <- "YifD"
    if (PFx >= PFx95h) Idc <- "Yes "
  } else if (Id == 1) { # type-1 changepoints
    if (PFx < PFx95l) {
      Idc <- "No  "
    } else if (PFx >= PFx95l & PFx < PFx95h) {
      Idc <- "?   "
    } else if (PFx >= PFx95h) Idc <- "Yes "
  }
  Idc
}
