#' Penalized maximal t (PMT) test with a reference series
#'
#' 以基准序列与参考序列的差值序列作检验，识别平均突变点（Wang 2008b）。
#' 与 [FindU()] 的区别是：检验对象为 `Y0 = base - reference`，且临界值表
#' 由 `readPTtable()` 生成。
#'
#' @param Bseries base series，文件路径或数据框。
#' @param Rseries reference series，文件路径或数据框。
#' @param output 输出文件前缀。
#' @param MissingValueCode 缺测值标记。
#' @param plev 置信水平，取 0.75、0.80、0.90、0.95、0.99、0.9999 之一。
#' @param Iadj 订正到哪一段：`0` 为最长段，`> Ns + 1` 视为最后一段。
#' @param Mq 分位数匹配的类别数；`0` 按最短段长度自动定，`1` 为均值订正。
#' @param GUI 为 `TRUE` 时只返回 `0`，不返回结果列表。
#' @param Ny4a 用于输出的年数参数。
#' @param is_plot 是否画图。
#' @param ... 保留参数。
#' @param .state 由 [Read.wRef] 生成的 `RHState`；提供后不再读取数据。
#'
#' @section 数据流与参数来源:
#' 输入数据只读取一次，其余参数均由 [RHState] 携带：
#'
#' - `Bseries`、`Rseries`、`MissingValueCode`、`plev` → `Read.wRef()`
#' - `N <- length(Y0)` ← `state@series$Y0`（两序列共同时段长度）
#' - `Y0`、`Ti`、`IY0`、`IY0flg`、`TiB`、`IYBflg`、`bdata`、`ori.bdata`、
#'   `owflg` ← `state@series`
#' - `Nt` ← `state@season`
#' - `phi`、`PTmax` ← `state@critical_values`，由
#'   `readPTtable(N, plev)` 按序列长度与置信水平生成
#'
#' 传入 `.state` 时跳过读取，而 `plev` 必须与建表时一致，由
#' `check_rh_state()` 校验。
#'
#' @return
#' A list with:
#'
#' - `fit` : `data.table`，合并参考序列后的拟合与订正结果
#' - `TP`  : `data.table`，突变点检验结果
#'
#' @seealso [Read.wRef()]、[FindUD.wRef()]、[StepSize.wRef()]
#' @export
FindU.wRef <- function(
    Bseries = NULL, Rseries = NULL, output, MissingValueCode = "-999.99",
    plev = 0.95, Iadj = 10000, Mq = 10, GUI = FALSE, Ny4a = 0,
    is_plot = TRUE, ..., .state = NULL) {
  mkdir(dirname(output))
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
  Nadj <- Ny4a * Nt # Y0 is Base-Ref on common period

  Pk0 <- Pk.PMT(N) # calculate penalty vector for target Y0 series
  oout <- PTK(Y0, Pk0, Nmin) # find 1st break point, then get 2 more, pick max(prob)
  # as official first break point
  I0 <- 0
  I2 <- oout$KPx
  I4 <- N

  oout1 <- PTKI0I2(Y0, I0, I2, Nmin)
  I1 <- oout1$Ic

  oout2 <- PTKI0I2(Y0, I2, I4, Nmin)
  I3 <- oout2$Ic

  oout3 <- PTKI0I2(Y0, I1, I3, Nmin)
  I02 <- oout3$Ic

  prob <- PTKIc(Y0, Pk0, I02)$prob
  if (I02 > 0 & prob >= plev) {
    Ns <- 1
    Ips <- c(I02, N)
    found <- search_pmt_candidates(Y0, Ips, Ns, Nmin, phi, PTmax, plev)
    Ips <- found$Ips
    Ns <- found$Ns
  } else {
    Ns <- 0
    Ips <- N
  }

  Ns.initial <- Ns
  oout <- Rphi(Y0, Ips, Ns, IY0flg)
  cor <- oout$cor
  corL <- oout$corl
  corU <- oout$corh

  # find first possible changepoint
  PTx95 <- getPTx95(cor, N - 1, phi, PTmax)
  PTx95L <- getPTx95(corL, N - 1, phi, PTmax)
  PTx95U <- getPTx95(corU, N - 1, phi, PTmax)

  I0 <- 1
  Ic <- I02
  otmp <- PTKIc(Y0, Pk0, I02)
  ofileIout <- paste(output, "_1Cs.txt", sep = "")
  file.create(ofileIout)
  
  # ofileMout <- paste(output, "_mCs.txt", sep = "")
  # ofileRout<-paste(output,"_Base_Ref.fitU",sep="")
  if (otmp$PTk < PTx95L) { # search for more changepoints
    fprintf("%s changepoints in Series %s %s \n",
      Ns, Bseries, paste("sample:(", sprintf("%1.0f", 1), " ", sprintf("%-4.4s", "YifD"),
        sprintf("%10.0f", 19000101), ")",
        sep = ""
      ), file = ofileIout
    )
    cat("PMT finds no Type-1 changepoints in the series!\n")
    Ns <- 0
    Ips <- N
    #   return()
  } else {
    # start searching for possible new breakpoints
    Ns <- 1
    Ips <- c(I02, N) # first break point settle down, start searching for more
    tt <- TRUE
    Niter <- 0
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
      # finished find new possible changepoints, stored in Ips0
      tt1 <- TRUE
      while (tt1) { # check Ips0 and insert new break points
        if (length(Ips0) == 0) {
          tt1 <- FALSE
        } else {
          PTx.mx <- (-9999)
          PTx95L.mx <- .001
          for (i in 1:length(Ips0)) { # search Ips0 series, find max prob
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
            PTx95L <- getPTx95(corL, I2 - I0 - 1, phi, PTmax)
            otmp <- PTKIc(Y0[(I0 + 1):I2], Pk1, Ic - I0)
            if (otmp$PTk < PTx95L) {
              Ips0[i] <- 0
            } #         if(otmp$prob<plev) Ips0[i]<-0
            else {
              if (otmp$PTk / PTx95L > PTx.mx / PTx95L.mx) {
                PTx.mx <- otmp$PTk
                Imx <- Ic
                inc <- i
                PTx95L.mx <- PTx95L
              }
            }
          }
          if (PTx.mx >= PTx95L.mx) {
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

    # check least significant changepoint
    tt0 <- TRUE
    tt <- TRUE
    while (tt0) {
      while (tt) {
        PTx.mn <- 9999.9
        PTx95L.mn <- 9999.9
        Imin <- 0
        if (Ns == 0) {
          tt <- FALSE
        } else {
          for (i in 1:Ns) {
            I1 <- if (i == 1) 0 else Ips[i - 1]
            I3 <- Ips[i + 1]
            Ic <- Ips[i]
            Nseg <- I3 - I1
            Pk0 <- Pk.PMT(Nseg)
            PTx95 <- getPTx95(corL, Nseg - 1, phi, PTmax)
            otmp <- PTKIc(Y0[(I1 + 1):I3], Pk0, Ic - I1)
            #     otmp<-PTKI0I2(W,I1,I3, Nmin)
            if (otmp$PTk / PTx95 < PTx.mn / PTx95L.mn) {
              PTx.mn <- otmp$PTk
              PTx95L.mn <- PTx95
              Imin <- i
            }
          }
          if (Imin > 0 & PTx.mn < PTx95L.mn) {
            Ns <- Ns - 1
            Ips <- Ips[-Imin]
          } else {
            tt <- FALSE
          }
        }
      }
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
        tt <- T
        corL <- corl
      } else {
        tt0 <- FALSE
      }
    }
  }
  # final output
  Iseg.adj <- adjustment_segment(Ips, Ns, Iadj)

  ofileSout <- paste(output, "_Ustat.txt", sep = "")
  file.create(ofileSout)
  fprintf("Input Base Series: %s \n", Bseries, file = ofileSout)
  fprintf("Input Ref Series: %s \n", Rseries, file = ofileSout, append = TRUE)

  d_TP <- list()
  if (Ns > 0) {
    for (i in 1:Ns) {
      I0 <- if (i == 1) 0 else Ips[i - 1]
      Ic <- Ips[i]
      I3 <- Ips[i + 1]
      Nseg <- I3 - I0
      PTx95 <- getPTx95(cor, Nseg, phi, PTmax)
      PTx95l <- getPTx95(corl, Nseg, phi, PTmax)
      PTx95h <- getPTx95(corh, Nseg, phi, PTmax)
      Pk0 <- Pk.PMT(Nseg)
      otmpW <- PTKIc(W[(I0 + 1):I3], Pk0, Ic - I0)
      otmpL <- PTKIc(WL[(I0 + 1):I3], Pk0, Ic - I0)
      otmpU <- PTKIc(WU[(I0 + 1):I3], Pk0, Ic - I0)
      probW <- otmpW$prob
      probL <- min(c(otmpL$prob, otmpU$prob))
      probU <- max(c(otmpL$prob, otmpU$prob))
      otmp <- PTKIc(Y0[(I0 + 1):I3], Pk0, Ic - I0)
      #     else if(Id==1) { # type-1 changepoints
      if (otmp$PTk < PTx95l) {
        Idc <- "No  "
      } else if (otmp$PTk >= PTx95l & otmp$PTk < PTx95h) {
        Idc <- "?   "
      } else if (otmp$PTk >= PTx95h) Idc <- "Yes "
      #     }
      fprintf("PMT : c=%s; (Time %s); Type= 1; p=%s(%s-%s); PTmax=%s; CV95=%s(%s-%s); Nseg=%s\n",
        sprintf("%4.0f", Ic), sprintf("%10.0f", IY0[Ic]), sprintf("%10.4f", probW),
        sprintf("%6.4f", probL), sprintf("%6.4f", probU), sprintf("%10.4f", otmp$PTk),
        sprintf("%10.4f", PTx95), sprintf("%10.4f", PTx95l), sprintf("%10.4f", PTx95h),
        sprintf("%4.0f", Nseg), file = ofileSout, append = TRUE
      )
      d_TP[[i]] <- data.table(
        kind = 1, Idc, date = IY0[Ic],
        probL, probU, plev, PFx = otmp$PTk, PFx95l = PTx95l, PFx95h = PTx95h
      )
    }
  }

  if (Ns == 0) {
    fprintf("%s changepoints in Series %s , %s %s \n",
      Ns, Bseries, Rseries, paste("sample:(", sprintf("%1.0f", 1), " ", sprintf("%-4.4s", "YifD"),
        sprintf("%10.0f", 19000101), ")",
        sep = ""
      ), file = ofileIout
    )
    cat("PMT finds no Type-1 changepoints in the series!\n")
  } else {
    fprintf("The adj-diff. autocor is: %.4f ( %.4f , %.4f p= %.4f )\n",
      cor, corl, corh, p.cor,
      file = ofileSout, append = TRUE)
    fprintf("%s changepoints in Series %s , %s \n", Ns, Bseries, Rseries,
      file = ofileIout, append = TRUE
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
    is_plot = is_plot, plot_fun = plot_FindU.ref, ofileSout = ofileSout
  )
  # write.table(file = paste0(output, "_U.dat"), odata, col.names = TRUE, row.names = F, na = MissingValueCode)

  if (GUI) {
    return(0)
  } else {
    # file.copy(from = ofileIout, to = ofileMout, overwrite = TRUE)
    # cat("FindU.wRef finished successfully...\n")
    list(fit = odata, TP = d_TP)
  }
}
