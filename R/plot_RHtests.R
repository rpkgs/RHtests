#' @import ggplot2
#' @export
plot_RHtests <- function(r, outfile = "RHtests-noref.pdf", ..., show=FALSE) {
  plot_output(r$data, outfile, show=show) 
  # write_fig(p, outfile, 10, 6, show=show)
  # if (is_plotly) plotly::ggplotly(p)
}

#' @export
plot_output <- function(data, outfile = NULL, ..., show=FALSE) {
  d2 <- melt(data, c("id", "date", "base"))
  # d2[value < -100, value := NA]
  p <- ggplot(d2, aes(date, value)) +
    geom_line(data = data, aes(y = base, color = NULL), color = alpha("black", 0.6)) +
    geom_line(aes(color = variable))
  
  if (!is.null(outfile)) {
    write_fig(p, outfile, 10, 6, show = show)
  } else p
}


plot_RHtests_multi <- function(obj, outfile = "RHtests.pdf", ..., show=FALSE) {
  dout <- map(obj$result, ~ .$data[, .(date = num2date(date), base, QM_adjusted)]) %>%
    melt_list("site")
  n <- length(obj$result)

  p <- ggplot(dout, aes(date, y = QM_adjusted - base)) +
    geom_line() +
    # geom_line(aes(date, QM_adjusted), color = "blue") +
    facet_wrap(~site, scales = "free", ncol = 2)
  write_fig(p, outfile, 10, 50 / 70 * n, show=show)
}
