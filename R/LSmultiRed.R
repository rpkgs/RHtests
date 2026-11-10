LSmultipleRed <- function(Y0, Ti, Ips, IY0flg, Nt) {
  Ns <- length(Ips) - 1
  N <- length(Y0)
  otmp <- LS_MultiBreaks(Y0, Ti, Ips)
  coef <- otmp$coef
  beta <- otmp$coef[2]
  resi <- otmp$resi
  otmp <- autocorlh(resi, IY0flg)
  cor <- otmp$cor
  corl <- otmp$corl
  corh <- otmp$corh
  resi <- resi + beta * Ti
  W1 <- resi / (1 - cor)
  W2 <- c(W1[1], (resi[2:N] - cor * resi[1:(N - 1)]) / (1 - cor))
  W <- c(1, IY0flg[1:(N - 1)]) * W2 + (!c(1, IY0flg[1:(N - 1)])) * W1

  otmp <- LS_OneBreak(W, Ti - mean(Ti), NA)
  beta <- otmp$coef[2]
  St0 <- sum((Ti - mean(Ti))^2)
  df <- (N - 2 - Ns - Nt)
  sigmaE2 <- otmp$SSE / df
  t.stat <- abs(beta) / sqrt(sigmaE2 / St0)
  p.tr <- pt(t.stat, df)
  betaL <- beta - qt(.975, df) * sqrt(sigmaE2 / St0)
  betaU <- beta + qt(.975, df) * sqrt(sigmaE2 / St0)
  itmp <- Y0 - beta * Ti
  mu <- rep(0, Ns + 1)
  meanhat <- mu
  for (i in 0:Ns) {
    I0 <- if (i == 0) 1 else Ips[i] + 1
    I2 <- Ips[i + 1]
    mu[i + 1] <- mean(itmp[I0:I2])
    meanhat[I0:I2] <- mu[i + 1] + beta * Ti[I0:I2]
    resi[I0:I2] <- Y0[I0:I2] - meanhat[I0:I2]
  }
  W1 <- resi
  W2 <- c(resi[1], resi[2:N] - cor * resi[1:(N - 1)])
  W3 <- c(resi[1], resi[2:N] - corl * resi[1:(N - 1)])
  W4 <- c(resi[1], resi[2:N] - corh * resi[1:(N - 1)])
  W <- c(1, IY0flg[1:(N - 1)]) * W2 + (!c(1, IY0flg[1:(N - 1)])) * W1
  WL <- c(1, IY0flg[1:(N - 1)]) * W3 + (!c(1, IY0flg[1:(N - 1)])) * W1
  WU <- c(1, IY0flg[1:(N - 1)]) * W4 + (!c(1, IY0flg[1:(N - 1)])) * W1
  for (i in 0:Ns) {
    I0 <- if (i == 0) 1 else Ips[i] + 1
    I2 <- Ips[i + 1]
    W[I0:I2] <- W[I0:I2] + mean(itmp[I0:I2]) + beta * Ti[I0:I2]
    WL[I0:I2] <- WL[I0:I2] + mean(itmp[I0:I2]) + beta * Ti[I0:I2]
    WU[I0:I2] <- WU[I0:I2] + mean(itmp[I0:I2]) + beta * Ti[I0:I2]
  }

  listk(W, WL, WU, coef, cor, corl, corh, resi, mu, meanhat,
    trend = beta, betaL, betaU, p.tr
  )
}

#' LSmultiRedCycle
#'
#' 估计含季节周期与多个均值突变的趋势项，序列长度、季节周期长度与突变点数
#' 都需显式传入，以避开全局变量。
#'
#' @param Y0 去季节后的序列。
#' @param Ti 时间序号（长度同 `Y0`）。
#' @param Ips 突变点位置，末位为序列长度。
#' @param Iseg.adj 订正目标段的序号。
#' @param IY0flg 相邻观测是否连续（0/1）。
#' @param seasonal_mean 各月日的季节均值，长度 `Nt`。
#' @param Imd 观测的月日（mmdd）。
#' @param Icy 季节周期内的月日序列，长度 `Nt`。
#' @param Nt 季节周期长度（年 1、月 12、日 365/366）。
#'
#' @return
#' - `resi`: residual
#' - `meanhat`: linear trend
#' - `seasonal_mean`: 拟合后的季节周期。
LSmultiRedCycle <- function(Y0, Ti, Ips, Iseg.adj, IY0flg, seasonal_mean, Imd, Icy, Nt) {
  N <- length(Y0)
  Ns <- length(Ips) - 1
  Niter <- 0
  tt <- TRUE
  current_cycle <- seasonal_mean
  while (tt) {
    tt <- FALSE
    Niter <- Niter + 1
    previous_cycle <- current_cycle

    otmp <- LSmultipleRed(Y0, Ti, Ips, IY0flg, Nt)
    trend <- otmp$trend
    betaL <- otmp$betaL
    betaU <- otmp$betaU
    resi <- otmp$resi
    cor <- otmp$cor
    corl <- otmp$corl
    corh <- otmp$corh
    p.tr <- otmp$p.tr
    meanhat <- otmp$meanhat
    mu <- otmp$mu
    W <- otmp$W
    WL <- otmp$WL
    WU <- otmp$WU

    if (Nt > 1) {
      itmp1 <- cbind(previous_cycle, Icy)
      itmp2 <- cbind(1:N, Imd)
      colnames(itmp2) <- c("idx", "Icy")
      itmp <- merge(itmp1, itmp2, by = "Icy")
      seasonal_values <- itmp[order(itmp[, "idx"]), "previous_cycle"]

      for (i in 1:(Ns + 1)) {
        I0 <- if (i == 1) 1 else Ips[i - 1] + 1
        I2 <- Ips[i]
        #       delta<-if(i==(Ns+1)) 0 else mu[i]-mu[Iseg.adj]
        delta <- mu[i] - mu[Iseg.adj]
        Y0[I0:I2] <- Y0[I0:I2] + seasonal_values[I0:I2] - delta
      }

      for (i in 1:Nt) current_cycle[i] <- mean(Y0[Imd == Icy[i]])
      cycle_sd <- sqrt(var(current_cycle))
      if (is.na(cycle_sd)) {
        tt <- FALSE
      } else {
        itmp1 <- cbind(current_cycle, Icy)
        itmp2 <- cbind(1:N, Imd)
        colnames(itmp2) <- c("idx", "Icy")
        itmp <- merge(itmp1, itmp2, by = "Icy")
        seasonal_values <- itmp[order(itmp[, "idx"]), "current_cycle"]

        for (i in 1:(Ns + 1)) {
          I0 <- if (i == 1) 1 else Ips[i - 1] + 1
          I2 <- Ips[i]
          #         delta<-if(i==(Ns+1)) 0 else mu[i]-mu[Iseg.adj]
          delta <- mu[i] - mu[Iseg.adj]
          Y0[I0:I2] <- Y0[I0:I2] - seasonal_values[I0:I2] + delta
        }

        max_change <- max(abs(current_cycle - previous_cycle))
        if (max_change > cycle_sd / 1000 & Niter < 20) tt <- TRUE
      }
    }
  }

  listk(trend, betaL, betaU,
    seasonal_mean = current_cycle, mu, cor, corl, corh, W, WL, WU, resi,
    Y0 = as.vector(Y0), meanhat = as.vector(meanhat), p.tr
  )
}
