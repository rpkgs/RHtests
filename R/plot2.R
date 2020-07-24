plot2 <- function(at, labels, ...) {
    plot(..., xaxt = "n")
    axis(side = 1, at = at, labels = labels)
}
