#' Test candidate changepoints with a reference series
#'
#' 对 [FindU.wRef()] 给出的候选突变点逐个检验：`kind = 0` 的点按概率判定，
#' `kind = 1` 的点按 PTmax 与临界值比较；不显著的点会被剔除。
#'
#' @inheritParams FindU.wRef
#' @param InCs 候选突变点，data.frame 或文件路径，列 `date`（YYYYMMDD）与 `kind`。
#' @param is_plot 是否绘图。
#'
#' @section 数据流与参数来源:
#' 同 [FindU.wRef()]：数据只由 `Read.wRef()` 读取一次，`Y0`、`IY0`、`IY0flg`
#' 取自 `state@series`，`Nt` 取自 `state@season`，`phi`/`PTmax` 取自
#' `state@critical_values`（由 `readPTtable(N, plev)` 生成）。
#'
#' @return
#' A list with:
#'
#' - `fit` : `data.table`，合并参考序列后的拟合与订正结果
#' - `TP`  : `data.table`，突变点检验结果，列 `kind`、`Idc`、`date`、`probL`、
#'   `probU`、`plev`、`PFx`、`PFx95l`、`PFx95h`
#'
#' @seealso [Read.wRef()]、[FindU.wRef()]、[StepSize.wRef()]
#' @export
FindUD.wRef <- function(
    Bseries = NULL, Rseries = NULL, InCs, output,
    MissingValueCode = "-999.99",
    plev = 0.95, Iadj = 10000, Mq = 10, GUI = F, Ny4a = 0,
    is_plot = TRUE, .state = NULL) {
  state <- if (is.null(.state)) {
    Read.wRef(Bseries, Rseries, MissingValueCode, plev)
  } else .state
  check_rh_state(state, reference = TRUE, plev = plev)
  c(Y0, IY0, IY0flg) %<-% state@series[c("Y0", "IY0", "IY0flg")]
  Nt <- state@season$Nt
  c(phi, PTmax) %<-% state@critical_values[c("phi", "PTmax")]

  Nmin <- 5
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
      message <- paste("FindUD.wRef: invalid changepoint indices:", paste(Ips, collapse = " "), "\n")
      stop(message, call. = FALSE)
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

  Ns.initial <- Ns
  oout <- Rphi(Y0, Ips, Ns, IY0flg)
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
        stop("PMT finds the series ", Bseries, " to be homogeneous!", call. = FALSE)
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
          PTx95L <- getPTx95(corL, I2 - I0, phi, PTmax)
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
    otmp <- Rphi(Y0, Ips, Ns, IY0flg)
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
      stop("PMT finds the series ", Bseries, " to be homogeneous!", call. = FALSE)
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
    PTx95 <- getPTx95(cor, Nseg, phi, PTmax)
    PTx95L <- getPTx95(corl, Nseg, phi, PTmax)
    PTx95U <- getPTx95(corh, Nseg, phi, PTmax)

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
  odata <- finalize_wref(
    series = state@series, Ips = Ips, Iseg.adj = Iseg.adj, Nt = Nt,
    output = output, Mq = Mq, Nadj = Nadj, Ny4a = Ny4a,
    is_plot = is_plot, plot_fun = plot_FindUD.ref, ofileSout = ofileSout
  )
  # write.table(odata, paste0(output, "_UD.dat"), col.names = TRUE, row.names = F, na = MissingValueCode)
  
  if (GUI) {
    return(0)
  } else {
    # file.copy(from = ofileIout, to = ofileMout, overwrite = TRUE)
    # cat("FindUD.wRef finished successfully...\n")
    list(fit = odata, TP = d_TP)
  }
}
