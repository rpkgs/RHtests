#' StepSize.wRef
#'
#' @example R/examples/run_ex02.R
#' @export
StepSize.wRef <- function(
  Bseries = NULL, Rseries = NULL, InCs, output,
  MissingValueCode = "-999.99",
  plev = 0.95, Iadj = 10000, Mq = 10, GUI = F, Ny4a = 0,
  is_plot = TRUE, ..., debug = TRUE
) {
  Read.wRef(Bseries, Rseries, MissingValueCode)

  if (Ny4a > 0 & Ny4a <= 5) Ny4a <- 5
  N <- length(Y0)
  Nadj <- Nt * Ny4a
  points <- stepsize_points(InCs, IY0)
  Ips <- points$Ips
  Ids <- points$Ids
  Ns <- points$Ns
  if (anyNA(Ips) || !identical(Ips, sort(Ips))) {
    ErrorMSG <<- paste(
      "StepSize.wRef: Ips read in from ", InCs, "error!\n",
      get("ErrorMSG", env = .GlobalEnv), "\n"
    )
    if (!GUI) cat(ErrorMSG)
    return(-1)
  }

  # if(Iadj==1|Iseg.longest==0) Iseg.adj<-Ns+1 else Iseg.adj<-Iseg.longest
  Iseg.adj <- stepsize_segment(Ips, Iadj)

  fit_diff <- Rphi(Y0, Ips, Ns)
  df <- N - 2 - Ns
  p.cor <- pt(abs(fit_diff$cor) * sqrt(df / (1 - fit_diff$cor^2)), df)
  ofileIout <- paste0(output, "_fCs.txt")
  ofileMout <- paste0(output, "_mCs.txt")
  ofilePdf <- paste0(output, "_F.pdf")
  ofileSout <- paste0(output, "_Fstat.txt")
  # ofileRout<-paste(output,"_Base_Ref.fitUDfinal",sep="")
  file.create(ofileSout)
  file.create(ofileIout)
  # file.create(ofileRout)

  if (debug) {
    cat(paste("Input Base Series:", Bseries, "\n"), file = ofileSout)
    cat(paste("Input Ref Series:", Rseries, "\n"), file = ofileSout, append = T)
    fprintf("The adj-diff. autocor is: %.4f ( %.4f , %.4f p= %.4f )\n",
      fit_diff$cor, fit_diff$corl, fit_diff$corh, p.cor, file = ofileSout, append = TRUE
    )
  }

  cat(paste(Ns, "changepoints in Series", Bseries, "\n"), file = ofileIout)
  d_TP <- stepsize_test_ref(Y0, IY0, Ips, Ids, fit_diff, plev, debug, ofileSout)

  # estimate delta from Y0 (Base-Ref)
  muDif <- rep(0, Ns + 1)
  segments <- stepsize_bounds(Ips)
  omuDif <- Y0
  oY0 <- Y0
  for (i in 1:(Ns + 1)) {
    I0 <- segments$start[i]
    I2 <- segments$end[i]
    muDif[i] <- mean(Y0[I0:I2])
    omuDif[I0:I2] <- muDif[i]
  }
  # write.table(cbind(IY0,round(oY0,4),round(omuDif,4),round(Wo,4)),
  #             file=ofileRout,col.names=F,row.names=F)

  # transfer Ips(Base-Ref) to Ips(Base)
  Ips0 <- Ips
  IY1 <- .make_date(bdata)
  IYM <- bdata[, 2] * 100 + bdata[, 3]
  inds <- sort(unique(IYM))
  rtmp <- cbind(1:length(IY1), IY1)

  Ips <- merge(IY0[Ips0], rtmp, by.x = 1, by.y = "IY1")[, 2]
  Ips[Ns + 1] <- length(IY1)
  segments <- stepsize_bounds(Ips)

  Ti <- TiB
  N <- length(TiB)
  IY0flg <- IYBflg

  cycle_base <- rmCycle(bdata)
  adjBase <- cycle_base$Base
  EB <- cycle_base$EB[match(IYM, inds)]

  if (Ns > 0) {
    for (i in 1:(Ns + 1)) {
      I0 <- segments$start[i]
      I2 <- segments$end[i]
      DeltaD <- muDif[i] - muDif[Iseg.adj]
      adjBase[I0:I2] <- adjBase[I0:I2] + EB[I0:I2] - DeltaD # adjBase is base series adj to Iseg.adj
    }
  }

  Ydata <- cbind(bdata[, 1:3], adjBase)
  cycle_adjusted <- rmCycle(Ydata) # adjusted and de-seasonalized Base series
  EB0 <- cycle_adjusted$EB
  EEBd <- mean(EB0)
  Aadj <- cycle_adjusted$Base # de-seasonalize Badj
  Ipd <- c(N)
  # dtmp<-LSmultipleRed(Aadj,Ti,Ipd)
  fit_adjusted <- LSmultiRedCycle(Aadj, Ti, Ipd, Iseg.adj, IY0flg, EB0, IYM)

  # de-seasonalized Base series
  Ipd <- length(cycle_base$Base)
  # dtmp   <- LSmultipleRed(tbase,Ti,Ipd)
  fit_base <- LSmultiRedCycle(cycle_base$Base, Ti, Ipd, Iseg.adj, IY0flg, cycle_base$EB, IYM)
  meanhat0 <- fit_base$meanhat
  Ehat0 <- mean(meanhat0)

  if (debug) {
    fprintf(
      "Ignore changepoints -> trend0 = %.6f ( %.6f , %.6f ) (p= %.4f ); cor= %.4f ( %.4f , %.4f )\n\n",
      fit_base$trend, fit_base$betaL, fit_base$betaU, fit_base$p.tr,
      fit_base$cor, fit_base$corl, fit_base$corh, file = ofileSout, append = TRUE
    )
    cat("Step-sizes estimated from difference series:\n",
      file = ofileSout, append = TRUE
    )

    if (Ns > 0) {
      cat(round(muDif[2:(Ns + 1)] - muDif[1:Ns], 4),
        file = ofileSout, append = TRUE, fill = 80
      )
    }
    fprintf(
      "\n after such adjustments, the base series trend= %.6f ( %.6f , %.6f ) (p= %.4f ); cor= %.3f ( %.3f , %.3f )\n\n",
      fit_adjusted$trend, fit_adjusted$betaL, fit_adjusted$betaU, fit_adjusted$p.tr,
      fit_adjusted$cor, fit_adjusted$corl, fit_adjusted$corh, file = ofileSout, append = TRUE
    )
  }

  EB <- cycle_base$EB[match(IYM, inds)]
  Base <- cycle_base$Base
  EPB <- mean(bdata[, "data"], na.rm = T)
  tt <- TRUE
  Niter <- 0
  while (tt) {
    Niter <- Niter + 1
    EB00 <- EB
    otmp <- LSmultipleRed(Base, Ti, Ips, IY0flg)
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
        I0 <- segments$start[i]
        I2 <- segments$end[i]
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
        I0 <- segments$start[i]
        I2 <- segments$end[i]
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
  Ehat <- mean(meanhat)
  meanhat0 <- meanhat0 - Ehat0 + Ehat
  QMout <- NULL
  if (Ns > 0) {
    #   Rb<-Base-otmp$trend*Ti+EB
    #   QMout<-QMadjGaussian(Rb,Ips,Mq,Iseg.adj,Nadj)
    QMout <- QMadjGaussian.wRef(bdata[, 4], bdata[, 4] - bdata[, 5], Ips, Mq, Iseg.adj, Nadj, Nt, Ny4a)
    B <- QMout$PA

    if (debug) {
      adjustment <- function(i) QMout$AdjM[i]
      # Fcat: frequency category boundaries; DP: differences in category means.
      stepsize_qm_log(QMout, segments, Iseg.adj, Ny4a, adjustment, ofileSout)
    }
  } else {
    B <- Base
  }
  # else B<-Base-otmp$trend*Ti+EB
  # B<-B+otmp$trend*Ti

  adj <- Base + EB
  adjB <- Base + EB
  meanhatD <- rep(0, N)
  if (Ns > 0) {
    for (i in 1:(Ns + 1)) {
      I1 <- segments$start[i]
      I2 <- segments$end[i]
      Delta <- sig[Iseg.adj + 1] - sig[i + 1]
      DeltaD <- muDif[Iseg.adj] - muDif[i]
      adj[I1:I2] <- adj[I1:I2] + Delta
      adjB[I1:I2] <- adjB[I1:I2] + DeltaD
      meanhatD[I1:I2] <- EEBd + fit_adjusted$mu[1] + fit_adjusted$trend * Ti[I1:I2] - DeltaD
    }
  } else {
    meanhatD <- EPB + fit_adjusted$mu[1] + fit_adjusted$trend * Ti
  }

  if (is_plot) {
    plot_StepSize.ref(
      output, Base, EB, EB1, B,
      oY0, omuDif, otmp, meanhatD,
      QMout, Mq, Ns, adj, adjB, Ips, Iseg.adj
    )
  }

  if (debug) {
    cat("Common trend TPR fit to the de-seasonalized Base series:\n",
      file = ofileSout, append = TRUE
    )
    fprintf("#steps=  %d ; trend= %.6f (p= %.4f ); cor= %.4f ( %.4f , %.4f )  p= %.4f \n",
      Ns, otrend, p.tro, corout[1], corout[2], corout[3], corout[4],
      file = ofileSout, append = TRUE
    )
  }

  if (Ns > 0) {
    stepsizes <- c()
    for (i in 1:Ns) {
      stepsizes[i] <- if (i == 1) sig[i + 2] else sig[i + 2] - sig[i + 1]
    }
    # stepsizes = c(sig[3], diff(sig)[-(1:2)])
    if (debug) cat(round(stepsizes, 4), file = ofileSout, append = TRUE, fill = 80)
  }

  ori.bdata <- as_tibble(ori.bdata)
  n <- nrow(ori.bdata)
  expand <- function(x) stepsize_expand(x, owflg, n)
  id <- rep(NA_real_, n)
  id[owflg] <- Ti
  odata <- tibble(
    id = id,
    date = .make_date(ori.bdata),
    base = as.double(ori.bdata$data),
    `base_trend+meanShift.ref` = expand(meanhatD),
    base_mean_adj.ref = expand(adjB),
    `base_trend+meanShift.base` = expand(otmp$meanhat + mean(EB1)),
    base_mean_adj.base = expand(adj),
    base_anomaly = expand(Base),
    fit_of_base_anomaly = expand(otmp$meanhat),
    `annualC+meanShifts` = expand(otmp$meanhat + EB),
    # if(QMout$Mq>1) # vital change here, kongdd, 2020-07-24
    QM_adjusted = if (Ns > 0) expand(B) else rep(NA_real_, n),
    fit_of_deseason_base = expand(meanhat0)
  )
  # odata[owflg,3] <- round(Base+EB,4)
  # odata[owflg,11]<-round(Ro,4)

  odata <- stepsize_leap(
    odata, ori.bdata, IY0[Ips0], Base - otmp$trend * Ti + EB,
    Ti, otmp$trend, QMout, Mq, Iseg.adj, Nadj, Nt,
    fitted_cols = c("base_trend+meanShift.ref", "base_trend+meanShift.base",
                    "fit_of_base_anomaly", "annualC+meanShifts", "fit_of_deseason_base"),
    adjusted_cols = c("base_mean_adj.ref", "base_mean_adj.base", "base_anomaly")
  )

  odata %<>% as.data.table()
  odata$date %<>% num2date()
  # write.table(file=paste0(output,"_F.dat"), odata, col.names=TRUE, row.names=F, na=MissingValueCode)

  d_TP %<>% do.call(rbind, .)
  d_TP$stepsize <- stepsizes
  d_TP[, 4:ncol(d_TP)] %<>% lapply(round, digits = 4)
  fwrite(d_TP, ofileIout, append = TRUE, col.names = TRUE)

  if (GUI) {
    return(0)
  } else {
    # file.copy(from=ofileIout,to=ofileMout,overwrite=TRUE)
    # cat("StepSize.wRef finished successfully...\n")
    list(data = odata, TP = d_TP)
  }
}

# Evaluate changepoints in the base-reference difference series.
stepsize_test_ref <- function(Y0, IY0, Ips, Ids, fit, plev, debug, ofileSout) {
  Ns <- length(Ips) - 1
  d_TP <- list()
  if (Ns > 0) {
    for (i in 1:Ns) {
      Ic <- Ips[i]
      Id <- Ids[i]
      I0 <- if (i == 1) 0 else Ips[i - 1]
      I3 <- Ips[i + 1] # from t-1:t+1
      Nseg <- I3 - I0
      PTx95 <- getPTx95(fit$cor, Nseg)
      PTx95L <- getPTx95(fit$corl, Nseg)
      PTx95U <- getPTx95(fit$corh, Nseg)

      Pk0 <- Pk.PMT(Nseg)
      otmp <- PTKIc(fit$W[(I0 + 1):I3], Pk0, Ic - I0)
      prob <- otmp$prob
      otmp <- PTKIc(fit$WL[(I0 + 1):I3], Pk0, Ic - I0)
      probL0 <- otmp$prob
      otmp <- PTKIc(fit$WU[(I0 + 1):I3], Pk0, Ic - I0)
      probU0 <- otmp$prob
      probL <- min(probL0, probU0)
      probU <- max(probL0, probU0)

      otmp <- PTKIc(Y0[(I0 + 1):I3], Pk0, Ic - I0)
      PTx0 <- otmp$PTk
      Idc <- stepsize_classify(Id, probL, probU, plev, PTx0, PTx95L, PTx95U)

      if (debug) {
        fprintf(
          paste0("PMT : c=%4d; (Time %10d); Type= %4d; p=%10.4f",
                 "(%10.4f-%10.4f); PTmax=%10.4f; CV95=%10.4f(%10.4f-%10.4f); Nseg=%4d\n"),
          Ic, IY0[Ic], Id, prob, probL, probU, PTx0, PTx95, PTx95L, PTx95U, Nseg,
          file = ofileSout, append = TRUE
        )
      }
      d_TP[[i]] <- data.table(
        kind = Id, Idc,
        date = IY0[Ic], Ic, Nseg, stepsize = NA,
        probL, probU, plev, PFx = PTx0, PFx95l = PTx95L, PFx95h = PTx95U, prob, PFx95 = PTx95
      )
    }
  }
  d_TP
}
