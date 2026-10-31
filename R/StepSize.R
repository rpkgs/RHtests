#' StepSize
#'
#' @param debug
#' - `true`  : `_Fstat.txt` will be shown
#' - `false` : no `_Fstat.txt`
#'
#' @example R/examples/run_ex01.R
#' @export
StepSize <- function(
  InSeries = NULL, InCs, output,
  MissingValueCode = "-999.99",
  GUI = FALSE, plev = 0.95, Iadj = 10000, Mq = 10, Ny4a = 0,
  is_plot = TRUE, ..., debug = TRUE
) {
  if (!is.null(InSeries)) data <- Read(InSeries, MissingValueCode) # data not used

  if (Ny4a > 0 & Ny4a <= 5) Ny4a <- 5
  N <- length(Y0)
  Nadj <- Ny4a * Nt

  points <- stepsize_points(InCs, IY0)
  Ips <- points$Ips
  Ids <- points$Ids
  Ns <- points$Ns
  if (anyNA(Ips) || !identical(Ips, sort(Ips))) {
    ErrorMSG <<- paste("StepSize: Ips read in from ", InCs, "error:")
    for (i in seq_len(Ns)) {
      ErrorMSG <<- paste(get("ErrorMSG", env = .GlobalEnv), Ips[i])
    }
    ErrorMSG <<- paste(get("ErrorMSG", env = .GlobalEnv), "\n\n")

    if (!GUI) cat(ErrorMSG)
    return(-1)
  }
  oout <- rmCycle(itable)
  Y1 <- oout$Base
  EB <- oout$EB
  if (length(EB) != length(Icy)) {
    ErrorMSG <<- paste(
      "Annual cycle length (from non-missing data) differ from original dataset",
      "\n", get("ErrorMSG", env = .GlobalEnv), "\n"
    )
    if (!GUI) cat(ErrorMSG)
    return(-1)
  }

  Iseg.adj <- stepsize_segment(Ips, Iadj)
  segments <- stepsize_bounds(Ips)

  Ip0 <- N
  fit_base <- LSmultiRedCycle(Y1, Ti, Ip0, 1, IY0flg, EB, Imd)
  meanhat0 <- fit_base$meanhat
  Ehat0 <- mean(meanhat0)

  otmp <- LSmultiRedCycle(Y1, Ti, Ips, Iseg.adj, IY0flg, EB, Imd)
  Y1 <- otmp$Y0
  cor <- otmp$cor
  corl <- otmp$corl
  corh <- otmp$corh
  pcor <- pt(abs(cor) * sqrt((N - 2) / (1 - cor^2)), N - 2)
  EB1 <- otmp$EB

  EBfull <- EB1[match(Imd, Icy)]
  EEB <- mean(EB1, na.rm = T)

  QMout <- NULL
  if (Ns > 0) {
    Rb <- Y1 - otmp$trend * Ti + EBfull
    QMout <- QMadjGaussian(Rb, Ips, Mq, Iseg.adj, Nadj, Nt)
    B <- QMout$PA
  } else {
    B <- Y1 - otmp$trend * Ti + EBfull
  }

  adj <- otmp$Y0 + EBfull
  B <- B + otmp$trend * Ti

  ofilePdf <- paste0(output, "_F.pdf")
  ofileSout <- paste0(output, "_Fstat.txt")
  ofileIout <- paste0(output, "_fCs.txt")
  ofileMout <- paste0(output, "_mCs.txt")
  file.create(ofileSout)
  file.create(ofileIout)
  file.create(ofilePdf)

  if (debug) {
    stepsize_log(
      ofileSout, plev, N, fit_base, Ns, QMout, Ny4a,
      Iseg.adj, Ips, otmp, cor, corl, corh, pcor
    )
  }

  stepsizes <- NULL
  if (Ns > 0) {
    for (i in 1:(Ns + 1)) {
      I1 <- segments$start[i]
      I2 <- segments$end[i]
      Delta <- otmp$mu[Iseg.adj] - otmp$mu[i]
      adj[I1:I2] <- adj[I1:I2] + Delta
      stepsize <- otmp$mu[i + 1] - otmp$mu[i]
      if (debug) {
        fprintf("%d %d stepsize= %.4f \n", Ips[i], IY0[Ips[i]], stepsize,
          file = ofileSout, append = TRUE
        )
      }
    }
    stepsizes <- otmp$mu[-1] - otmp$mu[-(Ns + 1)]
  }
  Ehat <- mean(otmp$meanhat)
  meanhat0 <- meanhat0 - Ehat0 + Ehat

  if (is_plot) {
    plot_stepsize(
      oout, ofilePdf, EBfull, EEB, B, QMout, Ms, Mq, Ns, adj,
      Ips, Iseg.adj, otmp
    )
  }

  odata <- stepsize_data(
    ori.itable, ooflg, otmp, EEB, adj, EBfull, meanhat0,
    Ns, QMout, B, Y1, Ti, Ips, IY0, Mq, Iseg.adj, Nadj, Nt
  )

  d_TP <- stepsize_test(
    Y1, Ti, Ips, Ids, IY0, otmp, plev, stepsizes, debug, ofileSout, ofileIout
  )

  odata %<>% as.data.table()
  odata$date %<>% num2date()
  # write.table(file = paste0(output, "_F.dat"), odata, na = MissingValueCode, col.names = TRUE, row.names = F)

  if (GUI) {
    return(0)
  } else {
    # file.copy(from=ofileIout, to=ofileMout,overwrite=TRUE)
    # cat("StepSize finished successfully...\n")
    list(data = as.data.table(odata), TP = d_TP) # fit = odata,
  }
}

# Write the fitted trend and adjustment summary.
stepsize_log <- function(ofileSout, plev, N, fit_base, Ns, QMout, Ny4a,
                         Iseg.adj, Ips, otmp, cor, corl, corh, pcor) {
  fprintf("The nominal level of confidence (1-alpha)= %.4f \n", plev, file = ofileSout)
  fprintf("Input data filename: N= %d \n", N, file = ofileSout, append = TRUE)
  fprintf(
    " Ignore changepoints -> trend0 = %.6f ( %.6f , %.6f ) (p= %.4f ); cor= %.4f ( %.4f , %.4f )\n",
    fit_base$trend, fit_base$betaL, fit_base$betaU, fit_base$p.tr,
    fit_base$cor, fit_base$corl, fit_base$corh, file = ofileSout, append = TRUE
  )
  cat("Common trend TPR fit to the seasonal-cycle-adjusted Base series:\n",
    file = ofileSout, append = T
  )

  if (Ns > 0) {
    segments <- stepsize_bounds(Ips)
    adjustment <- function(i) mean(QMout$W[segments$start[i]:segments$end[i]])
    # Fcat: frequency category boundaries; DP: differences in category means.
    stepsize_qm_log(QMout, segments, Iseg.adj, Ny4a, adjustment, ofileSout)
  }

  fprintf(
    "#steps=  %d ; trend= %.6f ( %.6f , %.6f ) (p= %.4f ); cor= %.4f ( %.4f , %.4f ) %.4f \n",
    Ns, otmp$trend, otmp$betaL, otmp$betaU, otmp$p.tr, cor, corl, corh,
    pcor, file = ofileSout, append = TRUE
  )
}

# Assemble fitted series and restore leap-day values.
stepsize_data <- function(ori.itable, ooflg, otmp, EEB, adj, EBfull, meanhat0,
                          Ns, QMout, B, Y1, Ti, Ips, IY0, Mq, Iseg.adj, Nadj, Nt) {
  ori.itable <- as_tibble(ori.itable)
  n <- nrow(ori.itable)
  expand <- function(x) stepsize_expand(x, ooflg, n)
  odata <- tibble(
    id = as.double(seq_len(n)),
    date = .make_date(ori.itable),
    base = as.double(ori.itable$data),
    `base_trend+meanShift` = expand(otmp$meanhat + EEB),
    base_mean_adj = expand(adj),
    base_anomaly = expand(otmp$Y0),
    fit_of_base_anomaly = expand(otmp$meanhat),
    `annualC+meanShifts` = expand(otmp$meanhat + EBfull),
    QM_adjusted = if (Ns > 0 && QMout$Mq > 1) expand(B) else rep(NA_real_, n),
    fit_of_deseason_base = expand(meanhat0)
  )
  # odata[ooflg,1] <- Ti
  # odata[ooflg,3] <- round(otmp$Y0+EBfull,4)
  # odata[ooflg,9] <- round(oR,4)

  odata <- stepsize_leap(
    odata, ori.itable, IY0[Ips], Y1 - otmp$trend * Ti + EBfull,
    Ti, otmp$trend, QMout, Mq, Iseg.adj, Nadj, Nt,
    fitted_cols = c(
      "base_trend+meanShift", "fit_of_base_anomaly",
      "annualC+meanShifts", "fit_of_deseason_base"
    ),
    adjusted_cols = c("base_mean_adj", "base_anomaly")
  )
  odata
}

# Evaluate each candidate changepoint and write its statistics.
stepsize_test <- function(Y1, Ti, Ips, Ids, IY0, fit, plev, stepsizes, debug,
                          ofileSout, ofileIout) {
  Ns <- length(Ips) - 1
  RW <- LSmultiple(fit$W, Ti, Ips)$resi
  RWL <- LSmultiple(fit$WL, Ti, Ips)$resi
  RWU <- LSmultiple(fit$WU, Ti, Ips)$resi

  d_TP <- list()
  if (Ns > 0) {
    for (i in 1:Ns) {
      I1 <- if (i == 1) 1 else Ips[i - 1] + 1
      I3 <- Ips[i + 1]
      Ic <- Ips[i]
      Id <- Ids[i]
      Nseg <- I3 - I1 + 1

      PFx95 <- getPFx95(fit$cor, Nseg)
      PFx95l <- getPFx95(fit$corl, Nseg)
      PFx95h <- getPFx95(fit$corh, Nseg)
      SSEf.Iseg <- sum(fit$resi[I1:I3]^2)
      Ips0 <- Ips[-i]
      otmp <- LSmultiple(Y1, Ti, Ips0)
      SSE0.Iseg <- sum(otmp$resi[I1:I3]^2)
      Fx <- (SSE0.Iseg - SSEf.Iseg) * (Nseg - 3) / SSEf.Iseg
      Pk0 <- Pk.PMFT(Nseg)
      PFx <- Fx * Pk0[Ic - I1 + 1]

      otmp <- LSmultiple(fit$W, Ti, Ips0)
      SSE0.Iseg <- sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg <- sum(RW[I1:I3]^2)
      Fx <- (SSE0.Iseg - SSEf.Iseg) * (Nseg - 3) / SSEf.Iseg
      if (Fx < 0) {
        PFx <- 0
        Fx <- 0
        prob <- 0
      } else {
        prob <- pf(Fx, 1, Nseg - 3)
      }

      otmp <- LSmultiple(fit$WL, Ti, Ips0)
      SSE0.Iseg <- sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg <- sum(RWL[I1:I3]^2)
      Fx <- (SSE0.Iseg - SSEf.Iseg) * (Nseg - 3) / SSEf.Iseg
      probL0 <- if (Fx < 0) 0 else pf(Fx, 1, Nseg - 3)

      otmp <- LSmultiple(fit$WU, Ti, Ips0)
      SSE0.Iseg <- sum(otmp$resi[I1:I3]^2)
      SSEf.Iseg <- sum(RWU[I1:I3]^2)
      Fx <- (SSE0.Iseg - SSEf.Iseg) * (Nseg - 3) / SSEf.Iseg
      probU0 <- if (Fx < 0) 0 else pf(Fx, 1, Nseg - 3)
      probL <- min(probL0, probU0)
      probU <- max(probL0, probU0)

      Idc <- stepsize_classify(Id, probL, probU, plev, PFx, PFx95l, PFx95h)

      # PMF: Regression
      if (debug) {
        fprintf(
          paste0("PMF : c=%4d; (Time %10d); Type= %4d; p=%10.4f",
                 "(%10.4f-%10.4f); PFmax=%10.4f; CV95=%10.4f(%10.4f-%10.4f); Nseg=%4d\n"),
          Ic, IY0[Ic], Id, prob, probL, probU, PFx, PFx95, PFx95l, PFx95h, Nseg,
          file = ofileSout, append = TRUE
        )
      }

      d_TP[[i]] <- data.table(
        kind = Id, Idc,
        date = IY0[Ic], Ic, Nseg, stepsize = stepsizes[i],
        probL, probU, plev, PFx, PFx95l, PFx95h, prob, PFx95
      )
    }

    if (Ns == 0) cat("PMF finds the series to be homogeneous!\n")
    fprintf("#  %d changepoints in Series \n", Ns, file = ofileIout)

    d_TP %<>% do.call(rbind, .)
    d_TP[, 6:13] %<>% lapply(round, digits = 4)
    fwrite(d_TP, ofileIout, sep = "\t", append = TRUE, col.names = TRUE)
    # fwrite(d_TP, ofileSout, sep = "\t", append = TRUE, col.names = TRUE)
  }
  d_TP
}
