#' @importFrom Ipaper is_empty rm_empty mkdir listk 
#' ok warn

range2 <- function(x, y, na.rm = TRUE) {
  c(min(x, y, na.rm = na.rm), max(x, y, na.rm = na.rm))
}


#' @export
num2date <- function(x) {
  if (is(x, "Date")) {
    return(x)
  }
  # else
  as.character(x) %>%
    gsub("-", "", .) %>%
    gsub("0000$", "0101", .) %>%
    gsub("00$", "01", .) %>%
    as.Date("%Y%m%d")
}

#' @export
date2num <- function(date) {
  as.character(date) %>%
    gsub("-", "", .) %>%
    as.numeric()
}

is.Date <- function(x) is(x, "Date")

plot2 <- function(at, labels, ...) {
  plot(..., xaxt = "n")
  axis(side = 1, at = at, labels = labels)
}
