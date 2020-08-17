library(R.matlab)
library(glue)
library(iterators)

l <- readMat("N:/Research/PML_V2/PML_V2/dat_Figure1.mat")
#
top5names <- unlist(l$top5names) %>% matrix(ncol = 5)
# write_list2xlsx(top5names, "a.xlsx")

val = array(1:length(top5names), dim = dim(top5names)) %>% t()

{
  n = 60
  cols = rep("white", n)
  col_grps = c("green", "skyblue", "yellow", "red")
  alpha_grps = c("Mean", "Std", "Seasonality")[-1]
  temp <- foreach(var = c("LAI", "P", "PET", "T"), i = icount()) %do% {
    I = grep(glue(" {var}"), top5names)
    cols[I] = col_grps[i]
  }
  I_sd = grep("Std", top5names)
  I_season = grep("Seasonality", top5names)
  cols[I_sd] %<>% alpha(0.5)
  cols[I_season] %<>% alpha(0.2)
}

environment(panel.levelplot.raster2) <- environment(panel.levelplot.raster)
panel.levelplot.raster2 <- function (x, y, z, subscripts, at = pretty(z), ...,
          zcol,
          col.regions = regions$col,
          alpha.regions = regions$alpha, interpolate = FALSE, identifier = "levelplot")
{
  if (length(subscripts) == 0)
    return()
  regions <- trellis.par.get("regions")
  x.is.factor <- is.factor(x)
  y.is.factor <- is.factor(y)
  x <- as.numeric(x)
  y <- as.numeric(y)
  z <- as.numeric(z)
  # zcol <- level.colors(z, at, col.regions, colors = TRUE)
  print(zcol)
  x <- x[subscripts]
  y <- y[subscripts]
  z <- z[subscripts]
  zcol <- zcol[subscripts]
  if (lattice:::hasGroupNumber())
    group <- list(...)$group.number
  else group <- 0
  if (x.is.factor) {
    ux <- seq(from = min(x, na.rm = TRUE), to = max(x, na.rm = TRUE))
    xwid <- 1L
  }
  else {
    ux <- sort(unique(x[!is.na(x)]))
    if (!isTRUE(all.equal(diff(range(diff(ux))), 0)))
      warning("'x' values are not equispaced; output may be wrong")
    xwid <- mean(diff(ux))
  }
  if (y.is.factor) {
    ux <- seq(from = min(y, na.rm = TRUE), to = max(y, na.rm = TRUE))
    ywid <- 1L
  }
  else {
    uy <- sort(unique(y[!is.na(y)]))
    if (!isTRUE(all.equal(diff(range(diff(uy))), 0)))
      warning("'y' values are not equispaced; output may be wrong")
    ywid <- mean(diff(uy))
  }
  ncolumns <- length(ux)
  nrows <- length(uy)
  xlow <- ux[1] - 0.5 * xwid
  xhigh <- ux[ncolumns] + 0.5 * xwid
  ylow <- uy[1] - 0.5 * ywid
  yhigh <- uy[nrows] + 0.5 * ywid
  zmat <- rep("transparent", ncolumns * nrows)
  idx <- match(x, ux)
  idy <- match(y, rev(uy))
  id <- idy + nrows * (idx - 1L)
  zmat[id] <- zcol
  dim(zmat) <- c(nrows, ncolumns)
  grid::grid.raster(as.raster(zmat), interpolate = interpolate, x = xlow,
              y = ylow, width = xhigh - xlow, height = yhigh - ylow,
              just = c("left", "bottom"), default.units = "native",
              name = trellis.grobname(paste(identifier, "raster",
                                            sep = "."), type = "panel", group = group))
}

{
  panel = function(x, y, z, subscripts, ...) {
    panel.levelplot.raster2(x, y, z, subscripts, zcol = cols)

    panel.abline(h = c(3, 6, 9)+0.5)
    panel.text(x, y, top5names)
  }
  p <- levelplot(val, panel = panel,
                ylab = "",
                xlab = expression(bold("Contribution rank")),
                colorkey = FALSE,
                scales = list(y = list(at = 1:12, cex = 0.95,
                                       labels = rep(c("Monthly", "Seasonal", "Yearly"), 4)),
                              cex = 0.95, fontface = 2),
                 aspect = 0.7)
  write_fig(p, "a.pdf", 10, 6)
}

