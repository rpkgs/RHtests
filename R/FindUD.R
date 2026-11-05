#' FindUD
#'
#' @inheritParams FindU
#'
#' @param InCs A data.frame of detected turning points returned from [FindU],
#' with the columns of `c("kind", "Idc", "date", "probL", "probU", "plev",
#' "PFx", "PFx95l", "PFx95h")`.
#'
#' @examples
#' f <- system.file("extdata/Example1.dat", package = "RHtests")
#' r1 <- FindU(f, output = file.path(tempdir(), "example01"))
#'
#' r2 <- FindUD(
#'   InSeries = f, InCs = r1$TP,
#'   output = file.path(tempdir(), "example01")
#' )
#' print(r2)
#'
#' plot_output(r2$fit)
#' @export
FindUD <- function(
  InSeries = NULL, InCs, output, MissingValueCode = "-999.99",
  GUI = FALSE, plev = 0.95, Iadj = 10000, Mq = 10, Ny4a = 0, is_plot = FALSE, .state = NULL
) {
  state <- if (is.null(.state)) {
    Read(InSeries, MissingValueCode, plev)
  } else {
    .state
  }
  check_rh_state(state, reference = FALSE, plev = plev)
  c(Y0, Ti, IY0, IY0flg, Imd, itable, ori.itable, ooflg) %<-%
    state@series[c("Y0", "Ti", "IY0", "IY0flg", "Imd", "itable", "ori.itable", "ooflg")]
  c(Icy, Nt) %<-% state@season[c("Icy", "Nt")]
  c(phi, PFmax) %<-% state@critical_values[c("phi", "PFmax")]

  Debug <- TRUE
  flog <- paste(output, ".log", sep = "")
  Nmin <- 10
  if (Ny4a > 0 && Ny4a <= 5) Ny4a <- 5

  N <- length(Y0)
  Nadj <- Ny4a * Nt

  # change InCs
  TP <- if (is.character(InCs)) fread(InCs) else InCs

  if (is.Date(TP$date)) TP$date %<>% date2num()
  Ns <- nrow(TP) # number of changing points
  if (is.null(Ns)) Ns <- 0

  Pk0 <- Pk.PMFT(N)
  oout <- RemoveSeasonal(Y0, Imd)
  Y1 <- oout$anomaly
  seasonal_mean <- oout$mu
  if (length(seasonal_mean) != length(Icy)) {
    stop("Annual cycle length differs from original dataset", call. = FALSE)
  }

  ofileIout <- paste(output, "_pCs.txt", sep = "")
  ofileMout <- paste(output, "_mCs.txt", sep = "")

  ofilePdf <- paste(output, "_UD.pdf", sep = "")
  ofileSout <- paste(output, "_UDstat.txt", sep = "")
  file.create(ofileIout)
  file.create(ofilePdf)
  file.create(ofileSout)
  fprintf("The nominal level of confidence (1-alpha)= %s \n", plev, file = ofileSout)
  fprintf("Input data filename: N= %s \n", N, file = ofileSout, append = TRUE)

  Ip0 <- N
  oout <- LSmultiRedCycle(Y1, Ti, Ip0, 1, IY0flg, seasonal_mean, Imd, Icy, Nt)
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
    oout <- PMFT(Y1, Ti, Pk0, Nmin)
    I0 <- 0
    I2 <- oout$KPx
    I4 <- N

    oout1 <- PMFxKxI0I2(Y1, Ti, I0, I2, Nmin)
    I1 <- oout1$Ic
    oout2 <- PMFxKxI0I2(Y1, Ti, I2, I4, Nmin)
    I3 <- oout2$Ic
    oout3 <- PMFxKxI0I2(Y1, Ti, I1, I3, Nmin)
    I2 <- oout3$Ic

    Ns <- 1
    indices <- c(I1, I2, I3)
    candidates <- pmf_candidates(indices, Y1, Ti, N, IY0flg)
    PFx <- vapply(candidates, `[[`, numeric(1), "PFx")
    best <- sort(PFx, decreasing = TRUE, index.return = TRUE)$ix[1]
    PFx.mx <- PFx[best]
    prob.mx <- candidates[[best]]$prob
    Imx <- indices[best]
    if (prob.mx < plev) {
      #     cat("PMF finds the series to be homogeneous!\n",file=ofileIout)
      fprintf("%s changepoints in Series \n", 0, file = ofileIout)
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
    Ips <- match_changepoint_dates(TP$date, IY0)
    Ids <- TP$kind

    if (sum(is.na(Ips)) > 0 || !identical(Ips, sort(Ips))) {
      message <- paste("FindUD: invalid changepoint indices:", paste(Ips, collapse = " "), "\n")
      stop(message, call. = FALSE)
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
        otmp <- PMFxKxI0I2(Y1, Ti, I0, I2, Nmin)
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
          for (i in seq_along(Ips0)) {
            Ips1 <- sort(c(Ips, Ips0[i]))
            ith <- match(Ips0[i], Ips1)
            otmp <- PMFxIseg(Y1, Ti, Ips1, ith, IY0flg, Nt, phi, PFmax)
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
    for (i in seq_along(Ips)) {
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
          otmp <- PMFxIseg(Y1, Ti, Ips, i, IY0flg, Nt, phi, PFmax)
          probL <- min(otmp$probL, otmp$probU)
          if (probL < probL.mn) {
            Iseg.mn <- i
            probL.mn <- probL
          }
        } # end if documented
      }
      if (Iseg.mn > 0 && probL.mn < plev) {
        Ips <- Ips[-Iseg.mn]
        Ids <- Ids[-Iseg.mn]
        Ns <- Ns - 1
        if (Ns > 0) tt <- TRUE
      }
    }
  }
  # all changepoints are significant, calculate stats and output

  Iseg.adj <- adjustment_segment(Ips, Ns, Iadj)

  otmp <- LSmultiRedCycle(Y1, Ti, Ips, Iseg.adj, IY0flg, seasonal_mean, Imd, Icy, Nt)
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
  EB1 <- otmp$seasonal_mean

  itmp1 <- cbind(EB1, Icy)
  itmp2 <- cbind(1:N, Imd)
  colnames(itmp2) <- c("idx", "Icy")

  itmp <- merge(itmp1, itmp2, by = "Icy")
  EBfull <- itmp[order(itmp[, "idx"]), "EB1"]
  EEB <- mean(EB1, na.rm = TRUE)

  if (Ns > 0) {
    Rb <- Y1 - otmp$trend * Ti + EBfull
    QMout <- QMadjGaussian(Rb, Ips, Mq, Iseg.adj, Nadj, Nt)
    B <- QMout$PA
    fprintf("Nseg_shortest = %s ; Mq =  %s ; Ny4a =  %s \n",
      QMout$Nseg.mn, QMout$Mq, Ny4a,
      file = ofileSout, append = TRUE
    )
    fprintf("\n Adjust to segment %s : from %s to %s \n",
      Iseg.adj, if (Iseg.adj == 1) 1 else Ips[Iseg.adj - 1] + 1,
      Ips[Iseg.adj],
      file = ofileSout, append = TRUE
    )
    #   cat("#Fcat, DP (CDF and Differnces in category mean)\n",file=ofileSout,
    #       append=TRUE)
    if (QMout$Mq > 1) {
      oline <- paste("#Fcat: frequency category boundaries\n",
        "#DP: Difference in the category means\n#",
        sep = ""
      )
      for (i in 1:Ns) {
        oline <- paste(
          oline, paste("Fcat.Seg", i, sep = ""),
          paste("DP.Seg", i, sep = "")
        )
      }

      oline <- paste(oline, "\n")
      cat(oline, file = ofileSout, append = TRUE)

      write.table(round(QMout$osmean, 4),
        file = ofileSout, append = TRUE, row.names = FALSE, col.names = FALSE
      )

      for (i in 1:(Ns + 1)) {
        I1 <- if (i == 1) 1 else Ips[i - 1] + 1
        I2 <- Ips[i]
        if (i != Iseg.adj) {
          fprintf("Seg. %s: mean of QM-adjustments =%.4f\n",
            i, mean(QMout$W[I1:I2]),
            file = ofileSout, append = TRUE
          )
        }
      }
    }
  } else {
    B <- Y1 - otmp$trend * Ti + EBfull
  }

  adj <- otmp$Y0 + EBfull
  # B<-B+EBfull+otmp$trend*Ti
  B <- B + otmp$trend * Ti

  fprintf(" Ignore changepoints -> trend0 = %s ( %s , %s ) (p= %s ); cor= %s ( %s , %s )\n",
    round(beta0, 6), round(betaL0, 6), round(betaU0, 6), round(p.tr0, 4), round(corD, 4),
    round(corDL, 4), round(corDU, 4),
    file = ofileSout, append = TRUE
  )
  cat("Common trend TPR fit to the seasonal-cycle-adjusted Base series:\n",
    file = ofileSout, append = TRUE
  )
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
      fprintf("%s %s stepsize= %.4f \n", Ips[i], IY0[Ips[i]],
        stepsize,
        file = ofileSout, append = TRUE
      )
    }
  }

  oR <- Y1 - otmp$meanhat
  oR[2:N] <- oR[2:N] - oR[1:(N - 1)] * cor
  Ehat <- mean(otmp$meanhat)
  meanhat0 <- meanhat0 - Ehat0 + Ehat

  if (is_plot) {
    plot_FindUD(
      oout, ofilePdf, EBfull, EEB, B, QMout, Mq, Ns, adj,
      Ips, Iseg.adj, otmp, ori.itable, ooflg
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
  if (sum(is.na(ori.itable[, 4]) == FALSE & Imd1 == 229) > 0) {
    if (Ns > 0) {
      tdata <- ori.itable[is.na(ori.itable[, 4]) == FALSE, ]
      IY1 <- tdata[, 1] * 10000 + tdata[, 2] * 100 + tdata[, 3]
      Ips.ymd <- IY0[Ips]
      Ips.1 <- rep(NA, Ns + 1)
      for (i in 1:Ns) Ips.1[i] <- c(seq_along(IY1))[IY1 == Ips.ymd[i]]
      Ips.1[Ns + 1] <- length(IY1)
      #     Ips.1<-c(seq_along(IY1))[Ips.ymd==IY1]
      Imd2 <- tdata[, 2] * 100 + tdata[, 3]
      Ids.leap <- c(seq_along(Imd2))[Imd2 == 229]
      Nl <- length(Ids.leap)
      Rb <- Y1 - otmp$trend * Ti + EBfull
      Rb1 <- tdata[, 4]
      Rb1[-Ids.leap] <- Rb
      Ti1 <- rep(NA, length(IY1))
      Ti1[-Ids.leap] <- Ti
      for (i in seq_along(Ids.leap)) {
        Rb1[Ids.leap[i]] <- tdata[Ids.leap[i], 4] + Rb1[Ids.leap[i] - 1] - tdata[Ids.leap[i] - 1, 4]
        Ti1[Ids.leap[i]] <- Ti1[Ids.leap[i] - 1]
      }
      if (QMout$Mq > 1) {
        B1 <- QMadjGaussian(Rb1, Ips.1, Mq, Iseg.adj, Nadj, Nt)$PA
        B1 <- B1 + otmp$trend * Ti1
        B1.leap <- B1[Ids.leap]
        odata[is.na(odata[, 3]) == FALSE & Imd1 == 229, 9] <- round(B1.leap, 4)
      }
    } else {
      odata[Imd1 == 229, 9] <- odata[Imd1 == 229, 3]
    }
    Ids.leapo <- c(1:dim(ori.itable)[1])[is.na(ori.itable[, 4]) == FALSE & Imd1 == 229]
    for (jth in seq_along(Ids.leapo)) {
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

  RW <- LS_MultiBreaks(W, Ti, Ips)$resi
  RWL <- LS_MultiBreaks(WL, Ti, Ips)$resi
  RWU <- LS_MultiBreaks(WU, Ti, Ips)$resi

  d_TP <- list()
  if (Ns > 0) {
    # d_TP = foreach(i = 1:Ns) %do% {
    for (i in 1:Ns) {
      I1 <- if (i == 1) 1 else Ips[i - 1] + 1
      I3 <- Ips[i + 1]
      Ic <- Ips[i]
      Id <- Ids[i]
      Nseg <- I3 - I1 + 1

      PFx95 <- getPFx95(cor, Nseg, phi, PFmax)
      PFx95l <- getPFx95(corl, Nseg, phi, PFmax)
      PFx95h <- getPFx95(corh, Nseg, phi, PFmax)
      SSEf.Iseg <- sum(Rf[I1:I3]^2)
      Ips0 <- Ips[-i]

      otmp <- LS_MultiBreaks(Y1, Ti, Ips0)
      SSE0.Iseg <- sum(otmp$resi[I1:I3]^2)
      Fx <- (SSE0.Iseg - SSEf.Iseg) * (Nseg - 3) / SSEf.Iseg

      Pk0 <- Pk.PMFT(Nseg)
      PFx <- Fx * Pk0[Ic - I1 + 1]

      otmp <- LS_MultiBreaks(W, Ti, Ips0)
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

      otmp <- LS_MultiBreaks(WL, Ti, Ips0)
      SSE0.Iseg <- sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg <- sum(RWL[I1:I3]^2)
      Fx <- (SSE0.Iseg - SSEf.Iseg) * (Nseg - 3) / SSEf.Iseg
      probL0 <- if (Fx < 0) 0 else pf(Fx, 1, Nseg - 3)

      otmp <- LS_MultiBreaks(WU, Ti, Ips0)
      SSE0.Iseg <- sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg <- sum(RWU[I1:I3]^2)
      Fx <- (SSE0.Iseg - SSEf.Iseg) * (Nseg - 3) / SSEf.Iseg
      probU0 <- if (Fx < 0) 0 else pf(Fx, 1, Nseg - 3)

      probL <- min(probL0, probU0)
      probU <- max(probL0, probU0)

      ## type_TP的判断标准
      # browser()
      Idc <- guess_sign_level(Id, plev, probL, probU, PFx, PFx95l, PFx95h)

      fprintf("%s %s%s (%s-%s)%s%s (%s-%s)\n",
        sprintf("%1.0f", as.numeric(Id)), sprintf("%-4.4s", Idc), sprintf("%10.0f", IY0[Ic]),
        sprintf("%6.4f", probL), sprintf("%6.4f", probU), sprintf("%6.3f", plev),
        sprintf("%10.4f", PFx), sprintf("%10.4f", PFx95l), sprintf("%10.4f", PFx95h),
        file = ofileIout, append = TRUE
      )
      fprintf("PMF : c=%s; (Time %s); Type= %s; p=%s(%s-%s); PFmax=%s; CV95=%s(%s-%s); Nseg=%s\n",
        sprintf("%4.0f", Ic), sprintf("%10.0f", IY0[Ic]), sprintf("%4.0f", as.numeric(Id)),
        sprintf("%10.4f", prob), sprintf("%10.4f", probL), sprintf("%10.4f", probU),
        sprintf("%10.4f", PFx), sprintf("%10.4f", PFx95), sprintf("%10.4f", PFx95l),
        sprintf("%10.4f", PFx95h), sprintf("%4.0f", Nseg),
        file = ofileSout, append = TRUE
      )
      d_TP[[i]] <- data.table(kind = Id, Idc, date = IY0[Ic], probL, probU, plev, PFx, PFx95l, PFx95h)
    }

    if (Ns == 0) cat("PMF finds the series to be homogeneous!\n")
    fprintf("#  %s changepoints in Series \n", Ns, file = ofileIout)

    if (!is_empty(d_TP)) {
      d_TP %<>% do.call(rbind, .)
      d_TP[, 4:9] %<>% lapply(round, digits = 4)
      fwrite(d_TP, ofileIout, append = TRUE, col.names = TRUE)
    }
  }

  odata %<>% as.data.table()
  odata$date %<>% num2date()
  # write.table(file = paste0(output, "_UD.dat"), odata, na = MissingValueCode, col.names = TRUE, row.names = FALSE)

  if (GUI) {
    return(0)
  } else {
    # file.copy(from = ofileIout, to = ofileMout, overwrite = TRUE)
    # cat("FindUD finished successfully...\n")
    list(fit = odata, TP = d_TP) # fit = odata,
  }
}

guess_sign_level <- function(Id, plev, probL, probU, PFx, PFx95l, PFx95h) {
  if (Id == 0) { # type-0 changepoints
    if (plev > probU) {
      Idc <- "No  "
    } else if (plev > probL && plev <= probU) {
      Idc <- "?   "
    } else if (plev <= probL) {
      Idc <- "YifD"
    }
    if (PFx >= PFx95h) Idc <- "Yes "
  } else if (Id == 1) { # type-1 changepoints
    if (PFx < PFx95l) {
      Idc <- "No  "
    } else if (PFx >= PFx95l && PFx < PFx95h) {
      Idc <- "?   "
    } else if (PFx >= PFx95h) {
      Idc <- "Yes "
    }
  }
  Idc
}
