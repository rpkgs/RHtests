QMadj.GaussianDLY <- function(
  InSeries, InCs, output, MissingValueCode, GUI = FALSE,
  Iadj = 10000, Mq = 10, Ny4a = 0
) {
  series <- ReadDLY.g(InSeries, MissingValueCode)
  Y0 <- series$Y0
  Ti <- series$Ti
  IY0 <- series$IY0
  IY0flg <- series$IY0flg
  Imd <- series$Imd
  Icy <- series$Icy
  Nt <- series$Nt
  itable <- series$itable
  ori.itable <- series$ori.itable
  olflg <- series$olflg
  N <- length(Y0)
  Nadj <- Ny4a * Nt
  itmp <- readLines(InCs)
  if (length(itmp) >= 2) {
    Ns <- length(itmp) - 1
    Ips <- c(rep(0, Ns), N)
    for (i in 1:Ns) { # using YYYYMMDD as index, searching for the largest
      # date less or equal to given YYYYMMDD
      ymdtmp <- as.numeric(substr(itmp[i + 1], 7, 16))
      if (ymdtmp == as.integer(ymdtmp / 100) * 100) ymdtmp <- ymdtmp + 15
      # set date as 15 if input break point is 00 for date
      it <- match(ymdtmp, IY0)
      if (!is.na(it)) {
        Ips[i] <- it
      } else {
        Ips[i] <- max(c(1:N)[IY0 <= ymdtmp])
      }
    }
    if (sum(is.na(Ips)) > 0 | !identical(Ips, sort(Ips))) {
      stop("QMadj.GaussianDLY: invalid changepoint indices: ", paste(Ips, collapse = " "), call. = FALSE)
    }
  } else {
    Ns <- 0
  }

  ofileSout <- paste(output, "_QMadjDLYstat.txt", sep = "")
  file.create(ofileSout)
  fprintf("Input data filename: ; N= %s \n", N, file = ofileSout)

  if (Ns > 0) {
    Nsegs <- Ips - c(0, Ips[1:Ns])
    Iseg.longest <- sort(Nsegs, index = TRUE, decreasing = TRUE)$ix[1]
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

  oout <- RemoveSeasonal(Y0, Imd)
  Y1 <- oout$anomaly
  EB <- oout$mu
  if (length(EB) != length(Icy)) {
    stop("Annual cycle length differs from original dataset", call. = FALSE)
  }

  otmp <- LSmultiRedCycle(Y1, Ti, Ips, Iseg.adj, IY0flg, EB, Imd, Icy, Nt)
  Y1 <- otmp$Y0
  cor <- otmp$cor
  corl <- otmp$corl
  corh <- otmp$corh
  pcor <- pt(abs(cor) * sqrt((N - 2) / (1 - cor^2)), N - 2)
  Rf <- otmp$resi
  W <- otmp$W
  L <- otmp$WL
  WU <- otmp$WU
  EB1 <- otmp$seasonal_mean

  itmp1 <- cbind(EB1, Icy)
  itmp2 <- cbind(1:N, Imd)
  colnames(itmp2) <- c("idx", "Icy")
  itmp <- merge(itmp1, itmp2, by = "Icy")
  EBfull <- itmp[order(itmp[, "idx"]), "EB1"]
  EEB <- mean(EB1, na.rm = TRUE)

  if (Ns > 0) {
    Rb <- Y1 - otmp$trend * Ti + EBfull
    #   QMout<-QMadjDLY(Rb,Ips,Mq,Iseg.adj)
    QMout <- QMadjGaussian(Rb, Ips, Mq, Iseg.adj, Nadj, Nt)
    B <- QMout$PA
  } else {
    B <- Y1 - otmp$trend * Ti + EBfull
  }

  adj <- otmp$Y0 + EBfull
  B <- B + otmp$trend * Ti

  if (Ns > 0) {
    fprintf("Nseg_shortest = %s ; Mq =  %s ; Ny4a =  %s \n",
      QMout$Nseg.mn, QMout$Mq, Ny4a,
      file = ofileSout, append = TRUE
    )
    fprintf("\n Adjust to segment %s : from %s to %s \n",
      Iseg.adj, if (Iseg.adj == 1) 1 else Ips[Iseg.adj - 1] + 1, Ips[Iseg.adj],
      file = ofileSout, append = TRUE
    )
    #   cat("#Fcat, DP (CDF and Differnces in category mean)\n",file=ofileSout,
    #       append=TRUE)
    if (QMout$Mq > 1) {
      oline <- paste("#Fcat: frequency category boundaries\n",
        "#DP: Difference in the category means\n#",
        sep = ""
      )
      for (i in 1:Ns) oline <- paste(oline, paste("Fcat.Seg", i, sep = ""), paste("DP.Seg", i, sep = ""))
      oline <- paste(oline, "\n")
      cat(oline, file = ofileSout, append = TRUE)

      write.table(round(QMout$osmean, 4),
        file = ofileSout, append = TRUE,
        row.names = FALSE, col.names = FALSE
      )
      for (i in 1:(Ns + 1)) {
        I1 <- if (i == 1) 1 else Ips[i - 1] + 1
        I2 <- Ips[i]
        if (i != Iseg.adj) {
          fprintf("Seg. %s: mean of QM-adjustments =%s\n",
            i, round(mean(QMout$W[I1:I2]), 4),
            file = ofileSout, append = TRUE
          )
        }
      }
    }
  }

  fprintf("#steps=  %s ; trend= %s ( %s , %s ) (p= %s ); cor= %s ( %s , %s ) %s \n",
    Ns, round(otmp$trend, 6), round(otmp$betaL, 6), round(otmp$betaU, 6), round(otmp$p.tr, 4),
    round(cor, 4), round(corl, 4), round(corh, 4), round(pcor, 4),
    file = ofileSout, append = TRUE
  )

  if (Ns > 0) {
    for (i in 1:(Ns + 1)) {
      I1 <- if (i == 1) 1 else Ips[i - 1] + 1
      I2 <- Ips[i]
      Delta <- otmp$mu[Iseg.adj] - otmp$mu[i]
      adj[I1:I2] <- adj[I1:I2] + Delta
      stepsize <- otmp$mu[i + 1] - otmp$mu[i]
      fprintf("%s %s stepsize= %s \n",
        Ips[i], IY0[Ips[i]], round(stepsize, 4),
        file = ofileSout, append = TRUE
      )
    }
  }

  # oR<-Y1-otmp$meanhat
  # oR[2:N]<-oR[2:N]-oR[1:(N-1)]*cor
  # Ehat<-mean(otmp$meanhat)

  if (Ns > 0) {
    odata <- matrix(NA, dim(ori.itable)[1], 8)
    odata[, 1] <- ori.itable[, 1]
    odata[, 2] <- ori.itable[, 2]
    odata[, 3] <- ori.itable[, 3]
    odata[olflg, 4] <- round(otmp$Y0 + EBfull, 4)
    if (QMout$Mq > 1) odata[olflg, 5] <- round(B, 4)
    odata[olflg, 6] <- round(adj, 4)
    if (QMout$Mq > 1) odata[olflg, 7] <- round(B - otmp$Y0 - EBfull, 4)
    odata[olflg, 8] <- round(adj - otmp$Y0 - EBfull, 4)
  } else {
    odata <- cbind(ori.itable[, c(1, 2, 3, 4, 4, 4)], 0, 0)
  }

  ofileAout <- paste(output, "_QMadjDLY.dat", sep = "")
  # ofileAout<-output # suggested by Lucie, user can choose whatever filename
  write.table(odata, file = ofileAout, na = MissingValueCode, col.names = FALSE, row.names = FALSE)

  ofilePdf <- paste(output, "_QMadjDLY.pdf", sep = "")
  pdf(file = ofilePdf, onefile = TRUE, paper = "letter")
  op <- par(no.readonly = TRUE) # the whole list of settable par's
  par(mfrow = c(2, 1), cex.main = .8, cex.lab = .8, cex.axis = .8, cex = .8)

  uyrs <- unique(floor(ori.itable[, 1] / 10)) * 10
  labels <- NULL
  ats <- NULL
  for (i in 1:length(uyrs)) {
    if (!is.na(match(uyrs[i], ori.itable[, 1]))) {
      labels <- c(labels, uyrs[i])
      ats <- c(ats, match(uyrs[i], ori.itable[, 1]))
    }
  }
  par(mar = c(3, 4, 3, 2) + .1)
  pdata <- rep(NA, dim(ori.itable)[1])
  pdata[olflg] <- otmp$Y0
  plot(1:dim(ori.itable)[1], pdata,
    type = "l", xlab = "", ylab = "",
    ylim = c(min(otmp$Y0, otmp$meanhat), max(otmp$Y0, otmp$meanhat)),
    xaxt = "n", col = "black", lwd = .3,
    main = "Base anomaly series and regression fit"
  )
  axis(side = 1, at = ats, labels = labels)
  pdata[olflg] <- otmp$meanhat
  lines(1:dim(ori.itable)[1], pdata, col = "red")

  pdata[olflg] <- otmp$Y0 + EBfull
  plot(1:dim(ori.itable)[1], pdata,
    type = "l", xlab = "", ylab = "",
    ylim = c(
      min(otmp$Y0 + EBfull, otmp$meanhat + EBfull),
      max(otmp$Y0 + EBfull, otmp$meanhat + EBfull)
    ),
    xaxt = "n", col = "black", lwd = .3,
    main = "Base series and regression fit"
  )
  axis(side = 1, at = ats, labels = labels)

  pdata[olflg] <- otmp$meanhat + EEB
  lines(1:dim(ori.itable)[1], pdata, col = "red")

  pdata[olflg] <- adj
  plot(1:dim(ori.itable)[1], pdata,
    type = "l", xlab = "", ylab = "",
    ylim = c(min(c(adj, B)), max(c(adj, B))),
    xaxt = "n", col = "black", lwd = .3,
    main = "Mean-adjusted base series"
  )
  axis(side = 1, at = ats, labels = labels)

  if (Ns > 0) {
    if (QMout$Mq > 1) {
      pdata[olflg] <- B
      plot(1:dim(ori.itable)[1], pdata,
        type = "l", xlab = "", ylab = "",
        ylim = c(min(c(adj, B)), max(c(adj, B))),
        xaxt = "n", col = "black", lwd = .3,
        main = "QM-adjusted base series"
      )
      axis(side = 1, at = ats, labels = labels)
    }
  }

  # QM 订正量分布图
  if (Ns > 0 && QMout$Mq > 1) plot_qm_adjustments(QMout, Ips, Iseg.adj, Ns)
  par(op)
  dev.off()
}

QMadj.GaussianDLY.wRef <- function(Bseries, Rseries, MissingValueCode, InCs, output, Iadj = 10000, Mq = 20, Ny4a = 5) {
  if (Ny4a > 0 & Ny4a <= 5) Ny4a <- 5

  state <- Read.wRef(Bseries, Rseries, MissingValueCode)
  series <- state@series
  Y0 <- series$Y0
  IY0 <- series$IY0
  IY0flg <- series$IY0flg
  TiB <- series$TiB
  IYBflg <- series$IYBflg
  bdata <- series$bdata
  ori.bdata <- series$ori.bdata
  owflg <- series$owflg
  Nt <- state@season$Nt
  N <- length(Y0)
  Nadj <- Nt * Ny4a
  itmp <- readLines(InCs)
  if (length(itmp) < 2) {
    #   stop("There is no input Ips")
    Ns <- 0
    Ips <- N
    Ids <- 0
  } else {
    Ns <- length(itmp) - 1
    Ips <- c(rep(0, Ns), N)
    Ids <- rep(0, Ns)
    for (i in 1:Ns) { # using YYYYMMDD as index, searching for the largest
      # date less or equal to given YYYYMMDD
      ymdtmp <- as.numeric(substr(itmp[i + 1], 7, 16))
      it <- match(ymdtmp, IY0)
      if (!is.na(it)) {
        Ips[i] <- it
      } else {
        Ips[i] <- max(c(1:N)[IY0 <= ymdtmp])
      }
      Ids[i] <- as.numeric(substr(itmp[i + 1], 1, 1))
    }
    if (sum(is.na(Ips)) > 0 | !identical(Ips, sort(Ips))) {
      stop("QMadj.GaussianDLY.wRef: invalid changepoint indices", call. = FALSE)
    }
  }
  if (Ns > 0) {
    Nsegs <- Ips - c(0, Ips[1:Ns])
    Iseg.longest <- sort(Nsegs, index = TRUE, decreasing = TRUE)$ix[1]
  } else {
    Iseg.longest <- 0
  }

  ofileSout <- paste(output, "_QMadjDLYwRefstat.txt", sep = "")
  file.create(ofileSout)
  fprintf("The base data filename is: %s ; N= %s \n", Bseries, N, file = ofileSout)
  fprintf("The ref data filename is: %s \n", Rseries, file = ofileSout, append = TRUE)

  # if(Iadj==1|Iseg.longest==0) Iseg.adj<-Ns+1 else Iseg.adj<-Iseg.longest
  if (Iadj > (Ns + 1) | Iseg.longest == 0) {
    Iseg.adj <- Ns + 1
  } else if (Iadj == 0) {
    Iseg.adj <- Iseg.longest
  } else {
    Iseg.adj <- Iadj
  }

  # estimate delta from Y0 (Base-Ref)
  otmp <- Rphi(Y0, Ips, Ns, IY0flg)
  cor <- otmp$cor
  muDif <- rep(0, Ns + 1)
  for (i in 1:(Ns + 1)) {
    I0 <- if (i == 1) 1 else Ips[i - 1] + 1
    I2 <- if (i > Ns) length(Y0) else Ips[i]
    muDif[i] <- mean(Y0[I0:I2])
  }

  # transfer Ips(Base-Ref) to Ips(Base)
  Ips0 <- Ips
  IY1 <- bdata[, 1] * 10000 + bdata[, 2] * 100 + bdata[, 3]
  IYM <- bdata[, 2] * 100 + bdata[, 3]
  inds <- sort(unique(IYM))
  rtmp <- cbind(1:length(IY1), IY1)
  Ips <- merge(IY0[Ips0], rtmp, by.x = 1, by.y = "IY1")[, 2]
  Ips[Ns + 1] <- length(IY1)

  Ti <- TiB
  N <- length(TiB)
  IY0flg <- IYBflg

  dtmp <- RemoveSeasonal(bdata[, 4], IYM)
  adjBase <- dtmp$anomaly
  EB <- rep(0, length(IY1))
  for (i in 1:length(IY1)) {
    EB[i] <- dtmp$mu[inds == IYM[i]]
  }

  if (Ns > 0) {
    for (i in 1:(Ns + 1)) {
      I0 <- if (i == 1) 1 else Ips[i - 1] + 1
      I2 <- if (i > Ns) N else Ips[i]
      DeltaD <- muDif[i] - muDif[Iseg.adj]
      adjBase[I0:I2] <- adjBase[I0:I2] + EB[I0:I2] - DeltaD # adjBase is base series adj to last seg
    }
  }
  EPBa <- mean(adjBase)

  dtmp <- RemoveSeasonal(adjBase, IYM) # adjusted and de-seasonalized Base series
  EB0 <- dtmp$mu
  EEBd <- mean(EB0)
  EB <- rep(0, length(IY1))
  for (i in 1:length(IY1)) {
    EB[i] <- dtmp$mu[inds == IYM[i]]
  }
  Aadj <- dtmp$anomaly # de-seasonalize Badj
  Ipd <- c(N)
  dtmp <- LSmultipleRed(Aadj, Ti, Ipd, IY0flg, Nt)
  # muD<-dtmp$mu[1]+EPBa-mean(bdata[,4])
  muD <- dtmp$mu[1]
  betaD <- dtmp$trend

  otmp <- RemoveSeasonal(bdata[, 4], IYM)
  EB <- rep(0, length(IY1))
  for (i in 1:length(IY1)) {
    EB[i] <- otmp$mu[inds == IYM[i]]
  }
  Base <- otmp$anomaly
  EPB <- mean(bdata[, "data"], na.rm = TRUE)
  tt <- TRUE
  Niter <- 0
  while (tt) {
    Niter <- Niter + 1
    EB00 <- EB
    otmp <- LSmultipleRed(Base, Ti, Ips, IY0flg, Nt)
    meanhat <- otmp$meanhat
    df <- (N - 2 - Nt - Ns)
    p.cor <- pt(abs(otmp$cor) * sqrt(df / (1 - otmp$cor^2)), df)
    otrend <- otmp$trend
    corout <- c(otmp$cor, otmp$corl, otmp$corh, p.cor)
    sig <- otmp$sig
    p.tro <- otmp$p.tr
    if (Ns == 0) {
      EB1 <- EB00
      tt <- FALSE
    } else {
      for (i in 1:(Ns + 1)) {
        I0 <- if (i == 1) 1 else Ips[i - 1] + 1
        I2 <- if (i > Ns) length(Base) else Ips[i]
        Delta <- sig[i + 1] - sig[Iseg.adj + 1]
        Base[I0:I2] <- Base[I0:I2] + EB00[I0:I2] - Delta
      }
      # re-estimate seasonal cycle using adjusted series:
      EB1 <- rep(0, length(inds))
      for (i in 1:length(inds)) {
        EB1[i] <- mean(Base[IYM == inds[i]], na.rm = TRUE)
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
        if (max(abs(EB0 - EB1)) < VEB / 1000 | Niter > 20) tt <- FALSE
      }
      EB0 <- EB1
    }
  }

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
    meanhatD <- EPB + muD + betaD * Ti
  }

  if (Ns > 0) {
    QMout <- QMadjGaussian.wRef(bdata[, 4], bdata[, 4] - bdata[, 5], Ips, Mq, Iseg.adj, Nadj, Nt, Ny4a)
    B <- QMout$PA
    fprintf("Nseg_shortest = %s ; Mq =  %s ; Ny4a =  %s \n",
      QMout$Nseg.mn, QMout$Mq, Ny4a,
      file = ofileSout, append = TRUE
    )
    fprintf("\n Adjust to segment %s : from %s to %s \n",
      Iseg.adj, if (Iseg.adj == 1) 1 else Ips[Iseg.adj - 1] + 1, Ips[Iseg.adj],
      file = ofileSout, append = TRUE
    )
    #   cat("#Fcat, DP (CDF and Differnces in category mean)\n",file=ofileSout,
    #       append=TRUE)
    if (QMout$Mq > 1) {
      oline <- paste("#Fcat: frequency category boundaries\n",
        "#DP: Difference in the category means\n#",
        sep = ""
      )
      for (i in 1:Ns) oline <- paste(oline, paste("Fcat.Seg", i, sep = ""), paste("DP.Seg", i, sep = ""))
      oline <- paste(oline, "\n")
      cat(oline, file = ofileSout, append = TRUE)

      write.table(round(QMout$osmean, 4),
        file = ofileSout, append = TRUE,
        row.names = FALSE, col.names = FALSE
      )
      for (i in 1:(Ns + 1)) {
        I1 <- if (i == 1) 1 else Ips[i - 1] + 1
        I2 <- Ips[i]
        if (i != Iseg.adj) {
          fprintf("Seg. %s: mean of QM-adjustments =%s\n",
            i, round(QMout$AdjM[i], 4),
            file = ofileSout, append = TRUE
          )
        }
      }
    }
  } else {
    B <- Base
  }
  # else B<-Base-otmp$trend*Ti+EB
  # B<-B+otmp$trend*Ti

  ofilePdf <- paste(output, "_QMadjDLYwRef.pdf", sep = "")
  pdf(file = ofilePdf, onefile = TRUE, paper = "letter")
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  par(mfrow = c(3, 1))
  par(mar = c(3, 4, 3, 2) + .1, cex.main = .8, cex.lab = .8, cex.axis = .8, cex = .8)

  uyrs <- unique(floor(ori.bdata[, 1] / 10)) * 10
  labels <- NULL
  ats <- NULL
  for (i in 1:length(uyrs)) {
    if (!is.na(match(uyrs[i], ori.bdata[, 1]))) {
      labels <- c(labels, uyrs[i])
      ats <- c(ats, match(uyrs[i], ori.bdata[, 1]))
    }
  }

  pdata <- rep(NA, dim(ori.bdata)[1])
  pdata[owflg] <- Base + EB
  plot(1:dim(ori.bdata)[1], pdata,
    type = "l", xlab = "", ylab = "",
    ylim = c(min(Base + EB), max(Base + EB)),
    xaxt = "n", col = "black", lwd = 0.5,
    main = "a. Base series and Mean-adjusted  base series"
  )
  axis(side = 1, at = ats, labels = labels)
  # pdata[owflg]<-otmp$meanhat+mean(EB1)
  # lines(1:dim(ori.bdata)[1],pdata,col="red")
  pdata[owflg] <- meanhatD
  lines(1:dim(ori.bdata)[1], pdata, col = "blue")
  pdata[owflg] <- adjB
  lines(1:dim(ori.bdata)[1], pdata, col = "red")

  if (Ns > 0) {
    if (QMout$Mq > 1) {
      pdata[owflg] <- B
      plot(1:dim(ori.bdata)[1], pdata,
        type = "l", xlab = "", ylab = "",
        ylim = c(min(B), max(B)),
        xaxt = "n", col = "black", lwd = 0.5,
        main = "b. QM-adjusted base series"
      )
      axis(side = 1, at = ats, labels = labels)
    }
  }

  # QM 订正量分布图
  if (Ns > 0 && QMout$Mq > 1) plot_qm_adjustments(QMout, Ips, Iseg.adj, Ns)
  par(op)
  dev.off()

  if (Ns > 0) {
    odata <- matrix(NA, dim(ori.bdata)[1], 8)
    odata[, 1] <- ori.bdata[, 1]
    odata[, 2] <- ori.bdata[, 2]
    odata[, 3] <- ori.bdata[, 3]
    odata[, 4] <- ori.bdata[, 4]
    if (QMout$Mq > 1) odata[owflg, 5] <- round(B, 4) # QMadjusted series
    odata[owflg, 6] <- round(adjB, 4) # Mean-adjusted series
    if (QMout$Mq > 1) odata[owflg, 7] <- round(B - ori.bdata[owflg, 4], 4)
    odata[owflg, 8] <- round(adjB - ori.bdata[owflg, 4], 4)
  } else {
    odata <- cbind(ori.itable[, c(1, 2, 3, 4, 4, 4)], 0, 0)
  }
  ofileAout <- paste(output, "_QMadjDLYwRef.dat", sep = "")
  write.table(odata, file = ofileAout, na = MissingValueCode, col.names = FALSE, row.names = FALSE)
  return(0)
}
