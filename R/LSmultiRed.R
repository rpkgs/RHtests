LSmultipleRed <- function(Y0, Ti, Ips) {
  Ns <- length(Ips) - 1
  N <- length(Y0)
  otmp <- LSmultiple(Y0, Ti, Ips)
  sig <- otmp$sig
  beta <- otmp$sig[2]
  resi <- otmp$resi
  otmp <- autocorlh(resi, IY0flg)
  cor <- otmp$cor
  corl <- otmp$corl
  corh <- otmp$corh
  resi <- resi + beta * Ti
  W1 <- resi / (1 - cor)
  W2 <- c(W1[1], (resi[2:N] - cor * resi[1:(N - 1)]) / (1 - cor))
  W <- c(1, IY0flg[1:(N - 1)]) * W2 + (!c(1, IY0flg[1:(N - 1)])) * W1

  otmp <- LSmatrix(W, Ti - mean(Ti), NA)
  beta <- otmp$sig[2]
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

  list(
    W       = W,
    WL      = WL,
    WU      = WU,
    sig     = sig,
    cor     = cor,
    corl    = corl,
    corh    = corh,
    resi    = resi,
    mu      = mu,
    meanhat = meanhat,
    trend   = beta,
    betaL   = betaL,
    betaU   = betaU,
    p.tr    = p.tr
  )
}

#' LSmultiRedCycle
#'
#' @return
#' - `resi`: residual
#' - `meanhat`: linear trend
LSmultiRedCycle <- function(Y0, Ti, Ips, Iseg.adj) {
  N <- length(Y0)
  Ns <- length(Ips) - 1
  Niter <- 0
  tt <- TRUE
  EB1 <- EB
  # global variables: Imd, Icy
  while (tt) {
    tt <- FALSE
    Niter <- Niter + 1
    EB0 <- EB1

    otmp <- LSmultipleRed(Y0, Ti, Ips)
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
      itmp1 <- cbind(EB0, Icy)
      itmp2 <- cbind(1:N, Imd)
      colnames(itmp2) <- c("idx", "Icy")
      itmp <- merge(itmp1, itmp2, by = "Icy")
      EBfull <- itmp[order(itmp[, "idx"]), "EB0"]

      for (i in 1:(Ns + 1)) {
        I0 <- if (i == 1) 1 else Ips[i - 1] + 1
        I2 <- Ips[i]
        #       delta<-if(i==(Ns+1)) 0 else mu[i]-mu[Iseg.adj]
        delta <- mu[i] - mu[Iseg.adj]
        Y0[I0:I2] <- Y0[I0:I2] + EBfull[I0:I2] - delta
      }

      for (i in 1:Nt) EB1[i] <- mean(Y0[Imd == Icy[i]])
      VEB <- sqrt(var(EB1))
      if (is.na(VEB)) {
        tt <- FALSE
      } else {
        itmp1 <- cbind(EB1, Icy)
        itmp2 <- cbind(1:N, Imd)
        colnames(itmp2) <- c("idx", "Icy")
        itmp <- merge(itmp1, itmp2, by = "Icy")
        EBfull <- itmp[order(itmp[, "idx"]), "EB1"]

        for (i in 1:(Ns + 1)) {
          I0 <- if (i == 1) 1 else Ips[i - 1] + 1
          I2 <- Ips[i]
          #         delta<-if(i==(Ns+1)) 0 else mu[i]-mu[Iseg.adj]
          delta <- mu[i] - mu[Iseg.adj]
          Y0[I0:I2] <- Y0[I0:I2] - EBfull[I0:I2] + delta
        }

        DEBmx <- max(abs(EB1 - EB0))
        if (DEBmx > VEB / 1000 & Niter < 20) tt <- TRUE
      }
    }
  }
  
  list(
    trend   = trend,
    betaL   = betaL,
    betaU   = betaU,
    EB      = EB1,
    mu      = mu,
    cor     = cor,
    corl    = corl,
    corh    = corh,
    W       = W,
    WL      = WL,
    WU      = WU,
    resi    = resi,
    Y0      = as.vector(Y0),
    meanhat = as.vector(meanhat),
    p.tr    = p.tr
  )
}

rmCycle <- function(idata) {
  tdata <- if (is.data.table(idata)) {
    cbind(idata[, month * 100 + day], idata$data)
  } else {
    cbind(idata[, 2] * 100 + idata[, 3], idata[, 4])
  }

  inds <- sort(unique(tdata[, 1]))
  nx <- length(inds)
  mu <- rep(0, nx)
  for (i in 1:nx) {
    ind <- tdata[, 1] == inds[i]
    mu[i] <- mean(tdata[ind, 2], na.rm = TRUE)
    tdata[ind, 2] <- tdata[ind, 2] - mu[i]
  }
  list(EB = mu, Base = tdata[, 2])
}

#' Multiple breakpoints piecewise Regression (same slope, different interception)
#'
#' @param Ips breakpoints
#'
#' @return
#' - `sig`   : coefficients
#' - `fitted`: fitted value
#' - `resi`  : residual of fitted value
LSmultiple <- function(Y, T, Ips, ...) {
  Nx <- length(Y)
  Ns <- length(Ips) - 1
  X <- as.matrix(Y)
  # D  <- rep(1,Nx)
  D <- cbind(1, T - mean(T))
  if (Ns >= 1) {
    for (i in 1:Ns) {
      tmp <- rep(0, Nx)
      tmp[(Ips[i] + 1):Ips[i + 1]] <- 1
      D <- cbind(D, tmp)
    }
  }
  lm_solve(D, X, ...)
}

#' single-breakpoint piecewise regression (same slope, different interception)
#'
#' @param Y The response vector
#' @param T The predictor vector
#'
#' @param Ic the position of change point, `[1, Ic]` and `[Ic + 1, end]` are
#' corresponding to the perid
#'
#' #' @return
#' - `sig`   : coefficients
#' - `fitted`: fitted value
#' - `resi`  : residual of fitted value
#'
#' - `SSE`   : sum of square error
#' @export
LSmatrix <- function(Y, T_anorm, Ic, ...) {
  Nx <- length(Y)
  X <- as.matrix(Y)
  # D  <- rep(1, Nx)
  # T_anorm = T - mean(T)
  D <- cbind(1, T_anorm)
  if (!is.na(Ic)) {
    D <- cbind(D, c(rep(0, Ic), rep(1, Nx - Ic)))
  }
  lm_solve(D, X, ...)
}

lm_solve <- function(D, X, only.SSE = FALSE, ...) {
  # microbenchmark::microbenchmark(
  #   sig  = solve(t(D) %*% D) %*% t(D) %*% X,
  #   sig3  = as.matrix(.lm.fit(D, X)$coefficients),
  #   sig2 = chol2inv(chol(t(D) %*% D)) %*% t(D) %*% X,
  #   times = 100
  # )
  l <- .lm.fit(D, X)
  # sig  <- as.matrix(l$coefficients)
  fitted <- D %*% l$coefficients
  # resi <- X-fitted
  SSE <- sum(l$residuals^2)
  if (!only.SSE) {
    list(
      sig = l$coefficients,
      fitted = as.vector(fitted),
      resi = l$residuals,
      SSE = SSE
    )
  } else {
    list(SSE = SSE)
  }
}
