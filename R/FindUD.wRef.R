#' @export
FindUD.wRef <- function(
    Bseries = NULL, Rseries = NULL, InCs, output,
    MissingValueCode = "-999.99",
    plev = 0.95, Iadj = 10000, Mq = 10, GUI = F, Ny4a = 0,
    is_plot = TRUE) {
  Read.wRef(Bseries, Rseries, MissingValueCode)

  back_Ti <- Ti
  back_IY0flg <- IY0flg
  on.exit({
    assign("Ti", back_Ti, envir = .GlobalEnv)
    assign("IY0flg", back_IY0flg, envir = .GlobalEnv)
  })

  Nmin <- 5
  assign("Nmin", Nmin, envir = .GlobalEnv)
  if (Ny4a > 0 & Ny4a <= 5) Ny4a <- 5

  N <- length(Y0)
  Nadj <- Nt * Ny4a

  # readin Ips
  TP <- if (is.character(InCs)) fread(InCs) else InCs
  if (is.Date(TP$date)) TP$date %<>% date2num()
  Ns <- nrow(TP) # number of changing points
  if (is.null(Ns)) Ns <- 0

  Pk0 <- Pk.PMT(N)
  ofileIout <- paste(output, "_pCs.txt", sep = "")
  file.create(ofileIout)
  if (Ns == 0) { # no input changepoints
    Pk0 <- Pk.PMT(N)
    oout <- PTK(Y0, Pk0, Nmin)
    I0 <- 0
    I2 <- oout$KPx
    I4 <- N
    I1 <- PTKI0I2(Y0, I0, I2, Nmin)$Ic
    I3 <- PTKI0I2(Y0, I2, I4, Nmin)$Ic
    I02 <- PTKI0I2(Y0, I1, I3, Nmin)$Ic
    prob <- PTKIc(Y0, Pk0, I02)$prob
    if (I02 > 0 & prob >= plev) {
      Ns <- 1
      Ips <- c(I02, N)
      Ids <- c(0, 1)
    } else {
      Ns <- 0
      Ips <- N
      Ids <- 0
    }
    Ns.ini <- Ns
    Ips.ini <- Ips
    Ids.ini <- Ids
  } else {
    Ips <- c(rep(0, Ns), N)
    Ids <- rep(0, Ns)
    for (i in 1:Ns) { # using YYYYMMDD as index, searching for the largest
      # date less or equal to given YYYYMMDD
      # ymdtmp<-as.numeric(substr(itmp[i+1],7,16))
      ymdtmp <- TP$date[i]
      it <- match(ymdtmp, IY0)
      if (!is.na(it)) {
        Ips[i] <- it
      } else {
        Ips[i] <- max(c(1:N)[IY0 <= ymdtmp])
      }
      # Ids[i]<-as.numeric(substr(itmp[i+1],1,1))
    }
    Ids <- TP$kind

    if (sum(is.na(Ips)) > 0 | !identical(Ips, sort(Ips))) {
      ErrorMSG <<- paste(
        "FindUD.wRef: Ips read in from ", InCs, "error!\n",
        get("ErrorMSG", env = .GlobalEnv), "\n"
      )
      if (!GUI) cat(ErrorMSG)
      return(-1)
    }
    Ns.ini <- Ns
    Ips.ini <- Ips
    Ids.ini <- Ids
  }

  Niter <- 0 # take cor==0, search for all possible changepoints
  tt <- TRUE
  while (tt) {
    Ns.old <- Ns
    Niter <- Niter + 1
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
          Ns1 <- Ns + 1
          Nseg <- I2 - I0
          PTx95 <- getPTx95(0, (Nseg - 1))
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

  Ns.initial <- Ns
  oout <- Rphi(Y0, Ips, Ns)
  cor <- oout$cor
  corL <- oout$corl
  corU <- oout$corh
  W <- oout$W
  WL <- oout$WL
  WU <- oout$WU

  Ns <- Ns.ini
  Ips <- Ips.ini
  Ids <- Ids.ini

  if (Ns == 0) { # search for 1st possible changepoint
    Pk0 <- Pk.PMT(N)
    prob1 <- PTKIc(W, Pk0, I02)$prob
    prob2 <- PTKIc(WL, Pk0, I02)$prob
    prob3 <- PTKIc(WU, Pk0, I02)$prob
    probU <- max(c(prob1, prob2, prob3))
    if (probU < plev) {
      #     cat("PMT finds the series to be homogeneous!\n",file=ofileIout)
      cat(paste(Ns, "changepoints in Series", Bseries, ",", Rseries, "\n"),
        file = ofileIout
      )
      if (!GUI) {
        cat("PMT finds the series to be homogeneous!\n")
        return()
      } else {
        ErrorMSG <<- paste("PMT finds the series", Bseries, "to be homogeneous!\n",
          get("ErrorMSG", env = .GlobalEnv), "\n",
          sep = ""
        )
        return(-1)
      }
    } else {
      Ns <- 1
      Ip <- c(I02, N)
      Id <- c(0, 1)
    }
  }

  Ips.i <- Ips
  Niter <- 0
  tt <- TRUE
  while (tt) { # condition on is there any more Bps to insert in Ips?
    Niter <- Niter + 1
    tt <- FALSE
    Ips0 <- NULL
    for (i in 1:(Ns + 1)) { # search for new break points
      I0 <- if (i == 1) 0 else Ips[i - 1]
      I2 <- Ips[i]
      otmp <- PTKI0I2(Y0, I0, I2, Nmin)
      if (otmp$prob > 0) Ips0 <- sort(c(Ips0, otmp$Ic))
    }
    tt1 <- TRUE
    while (tt1) { # check and insert new break points
      if (length(Ips0) == 0) {
        tt1 <- FALSE
      } else {
        probU.mx <- (-1)
        probL.mx <- (-1)
        PTx.mx <- (-9999)
        PTx95L.mx <- .001
        for (i in 1:length(Ips0)) { # search Ips1 series, find max prob
          Ips1 <- sort(c(Ips, Ips0[i]))
          Ic <- Ips0[i]
          id <- match(Ic, Ips1)
          if (id == length(Ips1)) {
            print(Ips1)
            print(id)
            stop("error in FindSteps")
          }
          I0 <- if (id == 1) 0 else Ips1[id - 1]
          I2 <- Ips1[id + 1]
          Ns1 <- Ns + 1
          Nseg <- I2 - I0
          Pk1 <- Pk.PMT(Nseg)
          PTx95L <- getPTx95(corL, I2 - I0)
          PTx <- PTKIc(Y0[(I0 + 1):I2], Pk1, Ic - I0)$PTk
          prob1 <- PTKIc(W[(I0 + 1):I2], Pk1, Ic - I0)$prob
          prob2 <- PTKIc(WL[(I0 + 1):I2], Pk1, Ic - I0)$prob
          prob3 <- PTKIc(WU[(I0 + 1):I2], Pk1, Ic - I0)$prob
          probU <- max(c(prob1, prob2, prob3))
          probL <- min(c(prob1, prob2, prob3))
          if (probU < plev) {
            Ips0[i] <- 0
          } else {
            if (PTx / PTx95L > PTx.mx / PTx95L.mx) {
              PTx.mx <- PTx
              PTx95L.mx <- PTx95L
              probL.mx <- probL
              probU.mx <- probL
              Imx <- Ic
              inc <- i
            }
          }
        }
        if (probU.mx >= plev) {
          Ips <- sort(c(Ips, Imx)) # insert new point into Ips
          Ns <- Ns + 1
          Ips0 <- Ips0[-inc] # exclude co-responding point in Ips0
          tt <- TRUE # continue search
        } else {
          tt1 <- FALSE
        } # finish inserting new points into Ips
        Ips0 <- Ips0[Ips0 != 0]
      }
    }
  } # finish search for new break points
  Ids0 <- rep(NA, length(Ips))
  for (i in 1:length(Ips)) {
    if (Ips[i] %in% Ips.i) {
      Ids0[i] <- Ids[Ips.i == Ips[i]]
    } else {
      Ids0[i] <- 0
    }
  }
  Ids <- Ids0

  tt <- TRUE
  tt0 <- TRUE
  while (tt0) {
    while (tt) { # delete changepoints which are not significant
      tt <- FALSE
      probL.mn <- 9999
      Iseg.mn <- 0
      for (i in 1:Ns) { # check all changepoints
        if (Ids[i] == 0) { # check those un-documented
          I0 <- if (i == 1) 0 else Ips[i - 1]
          I3 <- Ips[i + 1]
          Ic <- Ips[i]
          Nseg <- I3 - I0
          Pk0 <- Pk.PMT(Nseg)
          PTx <- PTKIc(Y0[(I0 + 1):I3], Pk0, Ic - I0)$PTk
          prob1 <- PTKIc(W[(I0 + 1):I3], Pk0, Ic - I0)$prob
          prob2 <- PTKIc(WL[(I0 + 1):I3], Pk0, Ic - I0)$prob
          prob3 <- PTKIc(WU[(I0 + 1):I3], Pk0, Ic - I0)$prob
          probU <- otmp$prob
          probL <- min(prob1, prob2, prob3)
          if (probL < probL.mn) {
            probL.mn <- probL
            Iseg.mn <- i
          }
        } # end if documented
      } # end of do-loop from 1 ~ (Ns+1)
      if (Iseg.mn > 0 & probL.mn < plev) {
        Ips <- Ips[-Iseg.mn]
        Ids <- Ids[-Iseg.mn]
        Ns <- Ns - 1
        if (Ns > 0) tt <- TRUE
      }
    } # end of do-while
    otmp <- Rphi(Y0, Ips, Ns)
    W <- otmp$W
    WL <- otmp$WL
    WU <- otmp$WU
    cor <- otmp$cor
    corl <- otmp$corl
    corh <- otmp$corh
    df <- (N - 2 - Ns)
    p.cor <- pt(abs(cor) * sqrt(df / (1 - cor^2)), df)
    if (Ns.initial > Ns) {
      Ns.initial <- Ns
    } else {
      tt0 <- FALSE
    }
  }

  # all changepoints in the list are significant, final estimates of step size
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

  ofileMout <- paste(output, "_mCs.txt", sep = "")
  ofileSout <- paste(output, "_UDstat.txt", sep = "")
  # ofileRout<-paste(output,"_Base_Ref.fitUD",sep="")
  file.create(ofileSout)
  # file.create(ofileRout)
  cat(paste("Input Base Series:", Bseries, "\n"), file = ofileSout)
  cat(paste("Input Ref Series:", Rseries, "\n"), file = ofileSout, append = T)
  if (Ns == 0) {
    cat(paste(Ns, "changepoints in Series", Bseries, ",", Rseries, "\n"),
      file = ofileIout
    )
    #   cat("PMT finds the series to be homogeneous!\n",
    #       file=ofileIout)
    if (!GUI) {
      cat("PMT finds the series to be homogeneous!\n")
      return()
    } else {
      ErrorMSG <<- paste("PMT finds the series", Bseries, "to be homogeneous!\n",
        get("ErrorMSG", env = .GlobalEnv), "\n",
        sep = ""
      )
      return(-1)
    }
  } else {
    cat(paste(
      "The adj-diff. autocor is:", round(cor, 4), "(", round(corl, 4),
      ",", round(corh, 4), "p=", round(p.cor, 4), ")\n"
    ), file = ofileSout, append = T)
    cat(paste(Ns, "changepoints in Series", Bseries, "\n"),
      file = ofileIout
    )
  }
  d_TP <- list()
  for (i in 1:Ns) {
    Ic <- Ips[i]
    Id <- Ids[i]
    I0 <- if (i == 1) 0 else Ips[i - 1]
    I3 <- Ips[i + 1]
    Nseg <- I3 - I0
    PTx95 <- getPTx95(cor, Nseg)
    PTx95L <- getPTx95(corl, Nseg)
    PTx95U <- getPTx95(corh, Nseg)

    Pk0 <- Pk.PMT(Nseg)
    otmp <- PTKIc(W[(I0 + 1):I3], Pk0, Ic - I0)
    prob <- otmp$prob
    otmp <- PTKIc(WL[(I0 + 1):I3], Pk0, Ic - I0)
    probL0 <- otmp$prob
    otmp <- PTKIc(WU[(I0 + 1):I3], Pk0, Ic - I0)
    probU0 <- otmp$prob
    probL <- min(probL0, probU0)
    probU <- max(probL0, probU0)

    otmp <- PTKIc(Y0[(I0 + 1):I3], Pk0, Ic - I0)
    PTx0 <- otmp$PTk
    if (Id == 0) { # type-0 changepoints
      if (probU < plev) {
        Idc <- "No  "
      } else if (probL < plev & probU >= plev) {
        Idc <- "?   "
      } else if (probL >= plev) Idc <- "YifD"
      if (PTx0 >= PTx95U) Idc <- "Yes "
    } else if (Id == 1) { # type-1 changepoints
      if (PTx0 < PTx95L) {
        Idc <- "No  "
      } else if (PTx0 >= PTx95L & PTx0 < PTx95U) {
        Idc <- "?   "
      } else if (PTx0 >= PTx95U) Idc <- "Yes "
    }

    cat(
      paste("PMT : c=", sprintf("%4.0f", Ic),
        "; (Time ", sprintf("%10.0f", IY0[Ic]),
        "); Type=", sprintf("%4.0f", as.numeric(Id)),
        "; p=", sprintf("%10.4f", prob),
        "(", sprintf("%10.4f", probL),
        "-", sprintf("%10.4f", probU),
        "); PTmax=", sprintf("%10.4f", PTx0),
        "; CV95=", sprintf("%10.4f", PTx95),
        "(", sprintf("%10.4f", PTx95L),
        "-", sprintf("%10.4f", PTx95U),
        "); Nseg=", sprintf("%4.0f", Nseg), "\n",
        sep = ""
      ),
      file = ofileSout, append = TRUE
    )
    d_TP[[i]] <- data.table(
      kind = Id, Idc, date = IY0[Ic],
      probL, probU, plev, PFx = otmp$PTk, PFx95l = PTx95L, PFx95h = PTx95U
    )
  }
  if (!is_empty(d_TP)) {
    d_TP %<>% do.call(rbind, .)
    d_TP[, 4:9] %<>% lapply(round, digits = 4)
    fwrite(d_TP, ofileIout, append = TRUE, col.names = TRUE)
  }

  # estimate delta and final output
  otmp <- Rphi(Y0, Ips, Ns)
  cor <- otmp$cor
  muDif <- rep(0, Ns + 1)
  Ro <- Y0
  Wo <- Y0
  omuDif <- Y0
  oY0 <- Y0
  for (i in 1:(Ns + 1)) {
    I0 <- if (i == 1) 1 else Ips[i - 1] + 1
    I2 <- if (i > Ns) length(Y0) else Ips[i]
    muDif[i] <- mean(Y0[I0:I2])
    omuDif[I0:I2] <- muDif[i]
    Ro[I0:I2] <- Y0[I0:I2] - muDif[i]
  }
  Wo[1] <- Ro[1]
  Wo[2:N] <- Ro[2:N] - cor * Ro[1:(N - 1)] * IY0flg[1:(N - 1)]
  # write.table(cbind(IY0,round(oY0,4),round(omuDif,4),round(Wo,4)),
  #             file=ofileRout,col.names=F,row.names=F)

  # transfer Ips(Base-Ref) to Ips(Base)
  Ips0 <- Ips
  IY1 <- bdata[, 1] * 10000 + bdata[, 2] * 100 + bdata[, 3]
  IYM <- bdata[, 2] * 100 + bdata[, 3]
  inds <- sort(unique(IYM))
  rtmp <- cbind(1:length(IY1), IY1)
  Ips <- merge(IY0[Ips0], rtmp, by.x = 1, by.y = "IY1", sort = T)[, 2]
  Ips[Ns + 1] <- length(IY1)

  Ti <- TiB
  assign("Ti", Ti, envir = .GlobalEnv)
  IY0flg <- IYBflg
  assign("IY0flg", IY0flg, envir = .GlobalEnv)
  N <- length(TiB)

  dtmp <- rmCycle(bdata)
  adjBase <- dtmp$Base
  EB <- rep(0, length(IY1))
  for (i in 1:length(IY1)) {
    EB[i] <- dtmp$EB[inds == IYM[i]]
  }
  for (i in 1:(Ns + 1)) {
    I0 <- if (i == 1) 1 else Ips[i - 1] + 1
    I2 <- if (i > Ns) N else Ips[i]
    DeltaD <- muDif[i] - muDif[Iseg.adj]
    adjBase[I0:I2] <- adjBase[I0:I2] + EB[I0:I2] - DeltaD # adjBase is base series adj to Iseg.adj
  }
  EPBa <- mean(adjBase)

  Ydata <- cbind(bdata[, 1:3], adjBase)
  dtmp <- rmCycle(Ydata) # adjusted and de-seasonalized Base series
  EB0 <- dtmp$EB
  EEBd <- mean(EB0)
  EB <- rep(0, length(IY1))
  for (i in 1:length(IY1)) {
    EB[i] <- dtmp$EB[inds == IYM[i]]
  }
  Aadj <- dtmp$Base # de-seasonalize Badj
  Ipd <- c(N)
  dtmp <- LSmultipleRed(Aadj, Ti, Ipd)
  # muD<-dtmp$mu[1]+EPBa-mean(bdata[,4])
  muD <- dtmp$mu[1]
  betaD <- dtmp$trend
  betaDL <- dtmp$betaL
  betaDU <- dtmp$betaU
  corD <- dtmp$cor
  corDL <- dtmp$corl
  corDU <- dtmp$corh
  p.trD <- dtmp$p.tr

  dtmp <- rmCycle(bdata) # de-seasonalized Base series
  tbase <- dtmp$Base
  Ipd <- length(tbase)
  dtmp <- LSmultipleRed(tbase, Ti, Ipd)
  beta0 <- dtmp$trend
  cor <- dtmp$cor
  corL <- dtmp$corl
  corU <- dtmp$corh
  meanhat0 <- dtmp$meanhat
  Ehat0 <- mean(meanhat0)

  cat(paste(
    "Ignore changepoints -> trend0 =", round(beta0, 6),
    "(", round(dtmp$betaL, 6), ",", round(dtmp$betaU, 6),
    ") (p=", round(dtmp$p.tr, 4),
    "); cor=", round(cor, 3), "(", round(corL, 3), ",",
    round(corU, 3), ")\n\n"
  ), file = ofileSout, append = TRUE)
  cat("Step-sizes estimated from difference series:\n",
    file = ofileSout, append = TRUE
  )
  cat(round(muDif[2:(Ns + 1)] - muDif[1:Ns], 4),
    file = ofileSout, append = TRUE, fill = 80
  )
  cat(paste(
    "\n after such adjustments, the base series trend=",
    round(betaD, 6), "(", round(betaDL, 6), ",", round(betaDU, 6),
    ") (p=", round(p.trD, 4), "); cor=",
    round(corD, 3), "(", round(corDL, 3),
    ",", round(corDU, 3), ")\n\n"
  ), file = ofileSout, append = TRUE)

  otmp <- rmCycle(bdata)
  EB <- rep(0, length(IY1))
  for (i in 1:length(IY1)) {
    EB[i] <- otmp$EB[inds == IYM[i]]
  }
  Base <- otmp$Base
  EPB <- mean(bdata[, "data"], na.rm = T)
  tt <- TRUE
  Niter <- 0
  while (tt) {
    Niter <- Niter + 1
    EB00 <- EB
    otmp <- LSmultipleRed(Base, Ti, Ips)
    p.cor <- pt(abs(otmp$cor) * sqrt((N - 2) / (1 - otmp$cor^2)), N - 2)
    meanhat <- otmp$meanhat
    corout <- c(otmp$cor, otmp$corl, otmp$corh, p.cor)
    sig <- otmp$sig
    p.tro <- otmp$p.tr
    for (i in 1:(Ns + 1)) {
      I0 <- if (i == 1) 1 else Ips[i - 1] + 1
      I2 <- if (i > Ns) length(Base) else Ips[i]
      Delta <- sig[i + 1] - sig[Iseg.adj + 1]
      Base[I0:I2] <- Base[I0:I2] + EB00[I0:I2] - Delta
    }
    # re-estimate seasonal cycle using adjusted series:
    EB1 <- rep(0, length(inds))
    for (i in 1:length(inds)) {
      EB1[i] <- mean(Base[IYM == inds[i]], na.rm = T)
    }
    for (i in 1:length(IY1)) {
      EB[i] <- EB1[inds == IYM[i]]
    }
    for (i in 1:(Ns + 1)) {
      I0 <- if (i == 1) 1 else Ips[i - 1] + 1
      I2 <- if (i > Ns) length(Base) else Ips[i]
      Delta <- sig[i + 1] - sig[Iseg.adj + 1]
      Base[I0:I2] <- Base[I0:I2] - EB[I0:I2] + Delta
    }
    VEB <- sqrt(var(EB1))
    if (is.na(VEB)) {
      tt <- FALSE
    } else {
      if (max(abs(EB0 - EB1)) < VEB / 1000 | Niter > 50) tt <- FALSE
    }
    EB0 <- EB1
  }
  Ehat <- mean(meanhat)
  meanhat0 <- meanhat0 - Ehat0 + Ehat
  Ro <- Base - meanhat
  Ro[2:N] <- Ro[2:N] - corout[1] * Ro[1:(N - 1)]

  if (Ns > 0) {
    #   Rb<-Base-otmp$trend*Ti+EB
    #   QMout<-QMadjGaussian(Rb,Ips,Mq,Iseg.adj,Nadj)
    QMout <- QMadjGaussian.wRef(bdata[, 4], bdata[, 4] - bdata[, 5], Ips, Mq, Iseg.adj, Nadj, Nt, Ny4a)
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
        "#DP: Difference in the category means\n#",
        sep = ""
      )
      for (i in 1:Ns) oline <- paste(oline, paste("Fcat.Seg", i, sep = ""), paste("DP.Seg", i, sep = ""))
      oline <- paste(oline, "\n")
      cat(oline, file = ofileSout, append = T)

      write.table(round(QMout$osmean, 4),
        file = ofileSout, append = T,
        row.names = F, col.names = F
      )
      for (i in 1:(Ns + 1)) {
        I1 <- if (i == 1) 1 else Ips[i - 1] + 1
        I2 <- Ips[i]
        if (i != Iseg.adj) {
          cat(paste("Seg. ", i, ": mean of QM-adjustments =", round(QMout$AdjM[i], 4),
            "\n",
            sep = ""
          ), file = ofileSout, append = T)
        }
      }
    }
  } else {
    B <- Base
  }
  # else B<-Base-otmp$trend*Ti+EB
  # B<-B+otmp$trend*Ti

  IY1 <- bdata[, 1] * 10000 + bdata[, 2] * 100 + bdata[, 3]
  adj <- Base + EB
  adjB <- Base + EB
  meanhatD <- rep(0, N)
  if (Ns > 0) {
    for (i in 1:(Ns + 1)) {
      I1 <- if (i == 1) 1 else Ips[i - 1] + 1
      I2 <- Ips[i]
      Delta <- sig[Iseg.adj + 1] - sig[i + 1]
      DeltaD <- muDif[Iseg.adj] - muDif[i]
      adj[I1:I2] <- adj[I1:I2] + Delta
      adjB[I1:I2] <- adjB[I1:I2] + DeltaD
      meanhatD[I1:I2] <- EEBd + muD + betaD * Ti[I1:I2] - DeltaD
    }
  } else {
    meanhatD <- EEBd + muD + betaD * Ti
  }

  if (is_plot) {
    plot_FindUD.ref(
      output, Base, EB, EB1, B,
      oY0, omuDif, otmp, meanhatD,
      QMout, Mq, Ns, adj, adjB, Ips, Iseg.adj
    )
  }

  cat("Common trend TPR fit to the de-seasonalized Base series:\n",
    file = ofileSout, append = TRUE
  )
  cat(
    paste(
      "#steps= ", Ns, "; trend=", round(sig[2], 6), "(p=",
      round(p.tro, 4), "); cor=",
      round(corout[1], 4), "(", round(corout[2], 4), ",", round(corout[3], 4),
      ")  p=", round(corout[4], 4), "\n"
    ),
    file = ofileSout, append = TRUE
  )
  oout <- NULL
  for (i in 1:Ns) {
    stepsize <- if (i == 1) sig[i + 2] else sig[i + 2] - sig[i + 1]
    oout <- c(oout, stepsize)
  }
  cat(round(oout, 4), file = ofileSout, append = TRUE, fill = 80)

  odata <- matrix(NA, dim(ori.bdata)[1], 12)
  odata[owflg, 1] <- Ti
  odata[, 2] <- ori.bdata[, 1] * 10000 + ori.bdata[, 2] * 100 + ori.bdata[, 3]
  # odata[owflg,3] <- round(Base+EB,4)
  odata[, 3] <- ori.bdata[, 4]
  odata[owflg, 4] <- round(meanhatD, 4)
  odata[owflg, 5] <- round(adjB, 4)
  odata[owflg, 6] <- round(otmp$meanhat + mean(EB1), 4)
  odata[owflg, 7] <- round(adj, 4)
  odata[owflg, 8] <- round(Base, 4)
  odata[owflg, 9] <- round(otmp$meanhat, 4)
  odata[owflg, 10] <- round(otmp$meanhat + EB, 4)
  # odata[owflg,11]<-round(Ro,4)
  if (Ns > 0) if (Mq > 1) odata[owflg, 11] <- round(B, 4)
  odata[owflg, 12] <- round(meanhat0, 4)

  Imd1 <- ori.bdata[, 2] * 100 + ori.bdata[, 3]
  if (sum(is.na(ori.bdata[, 4]) == F & Imd1 == 229) > 0) {
    if (Ns > 0) {
      tdata <- ori.bdata[is.na(ori.bdata[, 4]) == F, ]
      IY1 <- tdata[, 1] * 10000 + tdata[, 2] * 100 + tdata[, 3]
      Ips.ymd <- IY1[Ips]
      Ips.1 <- rep(NA, Ns + 1)
      for (i in 1:Ns) Ips.1[i] <- c(1:length(IY1))[IY1 == Ips.ymd[i]]
      Ips.1[Ns + 1] <- length(IY1)
      #     Ips.1<-c(1:length(IY1))[Ips.ymd==IY1]
      Imd2 <- tdata[, 2] * 100 + tdata[, 3]
      Ids.leap <- c(1:length(Imd2))[Imd2 == 229]
      Nl <- length(Ids.leap)
      Rb <- Base - otmp$trend * Ti + EB
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
        odata[is.na(odata[, 3]) == F & Imd1 == 229, 11] <- round(B1.leap, 4)
      }
    } else {
      odata[Imd1 == 229, 11] <- odata[Imd1 == 229, 3]
    }
    Ids.leapo <- c(1:dim(ori.bdata)[1])[is.na(ori.bdata[, 4]) == F & Imd1 == 229]
    for (jth in 1:length(Ids.leapo)) {
      kth <- Ids.leapo[jth]
      if (Ns > 0) {
        k1th <- if (odata[kth - 1, 2] %in% IY0[Ips]) (kth + 1) else (kth - 1)
      } else {
        k1th <- kth - 1
      }
      for (pth in c(4, 6, 9, 10, 12)) odata[kth, pth] <- odata[k1th, pth]
      for (pth in c(5, 7, 8)) {
        delta1 <- odata[k1th, 3] - odata[k1th, pth]
        odata[kth, pth] <- odata[kth, 3] - delta1
      }
    }
  }

  odata %<>% set_colnames(fitdata_varnames_ref) %>% data.table()
  odata$date %<>% num2date()
  # write.table(odata, paste0(output, "_UD.dat"), col.names = TRUE, row.names = F, na = MissingValueCode)
  
  if (GUI) {
    return(0)
  } else {
    # file.copy(from = ofileIout, to = ofileMout, overwrite = TRUE)
    # cat("FindUD.wRef finished successfully...\n")
    list(fit = odata, TP = d_TP)
  }
}
