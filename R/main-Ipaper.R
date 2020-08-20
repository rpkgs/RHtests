which.notnull <- function(x) {
    which(!sapply(x, is.null))
}
