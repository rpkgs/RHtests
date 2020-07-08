rm_empty <- function(x) {
    if (is.list(x)) {
        x[sapply(x, length) > 0]
    }
    else {
        x[!is.na(x)]
    }
}

range2 <- function(x, y, na.rm = TRUE) {
    c(min(x, y, na.rm = na.rm), max(x, y, na.rm = na.rm))
}
