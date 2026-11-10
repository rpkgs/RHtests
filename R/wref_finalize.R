# 由 FindU.wRef() 与 FindUD.wRef() 共用：已知差值序列的突变点位置后，
# 估计各段均值与基准序列趋势，把基准序列订正到 Iseg.adj 段，
# 追加写出统计信息，并组装对外输出表。
#
# 参数：
# - series   : RHState@series，提供 Y0（base - ref）、IY0、TiB、IYBflg、
#              bdata、ori.bdata、owflg
# - Ips      : 差值序列坐标下的突变点位置，末位为序列长度
# - Iseg.adj : 订正目标段的序号
# - plot_fun : plot_FindU.ref() 或 plot_FindUD.ref()
# 返回 data.table，date 已由 YYYYMMDD 转为 Date。
finalize_wref <- function(series, Ips, Iseg.adj, Nt, output, Mq, Nadj, Ny4a,
                          is_plot, plot_fun, ofileSout) {
  c(Y0, IY0, TiB, IYBflg, bdata, ori.bdata, owflg) %<-%
    series[c("Y0", "IY0", "TiB", "IYBflg", "bdata", "ori.bdata", "owflg")]
  Ns <- length(Ips) - 1L

  # 差值序列的分段均值：omuDif 把段均值展开到每个位置，供绘图使用。
  oY0 <- omuDif <- Y0
  muDif <- map_dbl(seq_len(Ns + 1L), \(i) {
    I0 <- if (i == 1) 1L else Ips[i - 1] + 1L
    mean(Y0[I0:(if (i > Ns) length(Y0) else Ips[i])])
  })
  for (i in seq_len(Ns + 1L)) {
    I0 <- if (i == 1) 1L else Ips[i - 1] + 1L
    I2 <- if (i > Ns) length(Y0) else Ips[i]
    omuDif[I0:I2] <- muDif[i]
  }

  # 把 Ips 由差值序列坐标转换到基准序列坐标。
  IY1 <- bdata[, 1] * 10000 + bdata[, 2] * 100 + bdata[, 3]
  IYM <- bdata[, 2] * 100 + bdata[, 3]
  inds <- sort(unique(IYM))
  rtmp <- cbind(seq_along(IY1), IY1)
  Ips <- merge(IY0[Ips], rtmp, by.x = 1, by.y = "IY1")[, 2]
  Ips[Ns + 1L] <- length(IY1)

  Ti <- TiB
  IY0flg <- IYBflg
  N <- length(IY1)

  # 先用差值序列的分段均值把基准序列订正到 Iseg.adj 段，得到订正后的趋势。
  cycle_base <- RemoveSeasonal(bdata$data, IYM)
  adjBase <- cycle_base$anomaly
  EB <- cycle_base$mu[match(IYM, inds)]
  for (i in seq_len(Ns + 1L)) {
    I0 <- if (i == 1) 1L else Ips[i - 1] + 1L
    I2 <- if (i > Ns) N else Ips[i]
    adjBase[I0:I2] <- adjBase[I0:I2] + EB[I0:I2] - (muDif[i] - muDif[Iseg.adj])
  }

  dtmp <- RemoveSeasonal(adjBase, IYM)
  EB0 <- dtmp$mu
  EEBd <- mean(EB0)
  EB <- EB0[match(IYM, inds)]
  dtmp <- LSmultipleRed(dtmp$anomaly, Ti, N, IY0flg, Nt)
  muD <- dtmp$mu[1]
  betaD <- dtmp$trend
  betaDL <- dtmp$betaL
  betaDU <- dtmp$betaU
  corD <- dtmp$cor
  corDL <- dtmp$corl
  corDU <- dtmp$corh
  p.trD <- dtmp$p.tr

  # 不考虑任何突变点时基准序列的趋势。
  tbase <- cycle_base$anomaly
  dtmp <- LSmultipleRed(tbase, Ti, length(tbase), IY0flg, Nt)
  meanhat0 <- dtmp$meanhat
  Ehat0 <- mean(meanhat0)

  fprintf("Ignore changepoints -> trend0 = %s ( %s , %s ) (p= %s ); cor= %s ( %s , %s )\n\n",
    round(dtmp$trend, 6), round(dtmp$betaL, 6), round(dtmp$betaU, 6), round(dtmp$p.tr, 4),
    round(dtmp$cor, 4), round(dtmp$corl, 4), round(dtmp$corh, 4), file = ofileSout, append = TRUE
  )
  cat("Step-sizes estimated from difference series:\n", file = ofileSout, append = TRUE)
  cat(round(diff(muDif), 4), file = ofileSout, append = TRUE, fill = 80)
  fprintf("\n after such adjustments, the base series trend= %s ( %s , %s ) (p= %s ); cor= %s ( %s , %s )\n\n",
    round(betaD, 6), round(betaDL, 6), round(betaDU, 6), round(p.trD, 4), round(corD, 3),
    round(corDL, 3), round(corDU, 3), file = ofileSout, append = TRUE
  )

  # 迭代：拟合多段趋势 → 订正 → 重估季节周期，直到季节周期稳定。
  EB <- cycle_base$mu[match(IYM, inds)]
  Base <- cycle_base$anomaly
  Niter <- 0L
  tt <- TRUE
  while (tt) {
    Niter <- Niter + 1L
    EB00 <- EB
    otmp <- LSmultipleRed(Base, Ti, Ips, IY0flg, Nt)
    meanhat <- otmp$meanhat
    df <- N - 2 - Ns
    p.cor <- pt(abs(otmp$cor) * sqrt(df / (1 - otmp$cor^2)), df)
    corout <- c(otmp$cor, otmp$corl, otmp$corh, p.cor)
    coef <- otmp$coef
    p.tro <- otmp$p.tr
    for (i in seq_len(Ns + 1L)) {
      I0 <- if (i == 1) 1L else Ips[i - 1] + 1L
      I2 <- if (i > Ns) N else Ips[i]
      Base[I0:I2] <- Base[I0:I2] + EB00[I0:I2] - (coef[i + 1] - coef[Iseg.adj + 1])
    }
    # 用订正后的序列重估季节周期，再把季节周期加回。
    EB1 <- map_dbl(inds, \(md) mean(Base[IYM == md], na.rm = TRUE))
    EB <- EB1[match(IYM, inds)]
    for (i in seq_len(Ns + 1L)) {
      I0 <- if (i == 1) 1L else Ips[i - 1] + 1L
      I2 <- if (i > Ns) N else Ips[i]
      Base[I0:I2] <- Base[I0:I2] - EB[I0:I2] + (coef[i + 1] - coef[Iseg.adj + 1])
    }
    VEB <- sqrt(var(EB1))
    tt <- !(is.na(VEB) || max(abs(EB0 - EB1)) < VEB / 1000 || Niter > 50)
    EB0 <- EB1
  }
  meanhat0 <- meanhat0 - Ehat0 + mean(meanhat)

  # 分位数匹配订正量（Mq = 1 即均值订正）。
  QMout <- NULL
  if (Ns > 0) {
    QMout <- QMadjGaussian.wRef(
      bdata[, 4], bdata[, 4] - bdata[, 5], Ips, Mq, Iseg.adj, Nadj, Nt, Ny4a
    )
    B <- QMout$PA
    fprintf("Nseg_shortest = %s ; Mq =  %s ; Ny4a =  %s \n",
      QMout$Nseg.mn, QMout$Mq, Ny4a, file = ofileSout, append = TRUE
    )
    fprintf("\n Adjust to segment %s : from %s to %s \n",
      Iseg.adj, if (Iseg.adj == 1) 1 else Ips[Iseg.adj - 1] + 1, Ips[Iseg.adj], file = ofileSout, append = TRUE
    )
    if (QMout$Mq > 1) {
      oline <- paste("#Fcat: frequency category boundaries\n",
        "#DP: Difference in the category means\n#",
        sep = ""
      )
      for (i in seq_len(Ns)) oline <- paste(oline, paste0("Fcat.Seg", i), paste0("DP.Seg", i))
      fprintf("%s \n", oline, file = ofileSout, append = TRUE)
      write.table(round(QMout$osmean, 4),
        file = ofileSout, append = TRUE, row.names = FALSE, col.names = FALSE
      )
      for (i in seq_len(Ns + 1L)) {
        if (i != Iseg.adj) {
          fprintf("Seg. %d: mean of QM-adjustments =%.4f\n", i, QMout$AdjM[i],
            file = ofileSout, append = TRUE
          )
        }
      }
    }
  } else {
    B <- Base
  }

  # 各段订正量：adj 用趋势段位移，adjB 用差值序列段位移。
  adj <- adjB <- Base + EB
  meanhatD <- EEBd + muD + betaD * Ti
  if (Ns > 0) {
    for (i in seq_len(Ns + 1L)) {
      ind <- if (i == 1) seq_len(Ips[1]) else seq(Ips[i - 1] + 1L, Ips[i])
      DeltaD <- muDif[Iseg.adj] - muDif[i]
      adj[ind] <- adj[ind] + (coef[Iseg.adj + 1] - coef[i + 1])
      adjB[ind] <- adjB[ind] + DeltaD
      meanhatD[ind] <- EEBd + muD + betaD * Ti[ind] - DeltaD
    }
  }

  if (is_plot) {
    # plot_FindUD.ref() 没有 bdata 形参，多余参数经它的 ... 忽略。
    plot_fun(
      output = output, Base = Base, EB = EB, EB1 = EB1, B = B,
      oY0 = oY0, omuDif = omuDif, otmp = otmp, meanhatD = meanhatD,
      QMout = QMout, Mq = Mq, Ns = Ns, adj = adj, adjB = adjB,
      Ips = Ips, Iseg.adj = Iseg.adj, IY0 = IY0, bdata = bdata,
      ori.bdata = ori.bdata, owflg = owflg
    )
  }

  cat("Common trend TPR fit to the de-seasonalized Base series:\n",
    file = ofileSout, append = TRUE
  )
  fprintf("#steps=  %s ; trend= %s (p= %s ); cor= %s ( %s , %s )  p= %s \n",
    Ns, round(coef[2], 6), round(p.tro, 4), round(corout[1], 4), round(corout[2], 4),
    round(corout[3], 4), round(corout[4], 4), file = ofileSout, append = TRUE
  )
  # 段间位移：首段取 coef[3]，其后取相邻系数之差。
  if (Ns > 0) {
    cat(round(c(coef[3], diff(coef[3:(Ns + 2)])), 4),
      file = ofileSout, append = TRUE, fill = 80
    )
  }

  # 输出表：ori.bdata 的行数，12 列（见 GLOBAL.R 的 fitdata_varnames_ref）。
  odata <- matrix(NA, nrow(ori.bdata), 12) %>% set_colnames(fitdata_varnames_ref)
  odata[owflg, 1] <- Ti
  odata[, 2] <- ori.bdata[, 1] * 10000 + ori.bdata[, 2] * 100 + ori.bdata[, 3]
  odata[, 3] <- ori.bdata[, 4]
  odata[owflg, 4] <- round(meanhatD, 4)
  odata[owflg, 5] <- round(adjB, 4)
  odata[owflg, 6] <- round(otmp$meanhat + mean(EB1), 4)
  odata[owflg, 7] <- round(adj, 4)
  odata[owflg, 8] <- round(Base, 4)
  odata[owflg, 9] <- round(otmp$meanhat, 4)
  odata[owflg, 10] <- round(otmp$meanhat + EB, 4)
  if (Ns > 0 && QMout$Mq > 1) odata[owflg, 11] <- round(B, 4) # QM
  odata[owflg, 12] <- round(meanhat0, 4)

  # 闰日 2-29 未参与拟合，其订正量按相邻观测的差值补齐。
  Imd1 <- ori.bdata[, 2] * 100 + ori.bdata[, 3]
  if (any(!is.na(ori.bdata[, 4]) & Imd1 == 229)) {
    if (Ns > 0) {
      tdata <- ori.bdata[!is.na(ori.bdata[, 4]), ]
      IY1 <- tdata[, 1] * 10000 + tdata[, 2] * 100 + tdata[, 3]
      Ids.leap <- which(tdata[, 2] * 100 + tdata[, 3] == 229)
      Rb <- Base - otmp$trend * Ti + EB
      Rb1 <- Ti1 <- rep(NA_real_, length(IY1))
      Rb1[-Ids.leap] <- Rb
      Ti1[-Ids.leap] <- Ti
      for (k in Ids.leap) {
        Rb1[k] <- tdata[k, 4] + Rb1[k - 1] - tdata[k - 1, 4]
        Ti1[k] <- Ti1[k - 1]
      }
      if (QMout$Mq > 1) {
        B1 <- QMadjGaussian(Rb1, Ips, Mq, Iseg.adj, Nadj, Nt)$PA + otmp$trend * Ti1
        odata[!is.na(odata[, 3]) & Imd1 == 229, 11] <- round(B1[Ids.leap], 4)
      }
    } else {
      odata[Imd1 == 229, 11] <- odata[Imd1 == 229, 3]
    }
    for (kth in which(!is.na(ori.bdata[, 4]) & Imd1 == 229)) {
      k1th <- if (Ns > 0 && odata[kth - 1, 2] %in% IY0[Ips]) kth + 1L else kth - 1L
      odata[kth, c(4, 6, 9, 10, 12)] <- odata[k1th, c(4, 6, 9, 10, 12)]
      for (pth in c(5, 7, 8)) {
        odata[kth, pth] <- odata[kth, 3] - (odata[k1th, 3] - odata[k1th, pth])
      }
    }
  }

  dt_num2date(odata)
}
