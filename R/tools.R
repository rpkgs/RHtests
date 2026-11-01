#' @importFrom Ipaper is_empty rm_empty mkdir ok warn

range2 <- function(x, y, na.rm = TRUE) {
  c(min(x, y, na.rm = na.rm), max(x, y, na.rm = na.rm))
}


#' 日期与 YYYYMMDD 数值互转
#'
#' 序列日期在内部统一为 YYYYMMDD 数值（月序列日为 `00`）。
#'
#' @param x,date 日期向量，`Date` 或可被 `as.Date()` 解析的字符串。
#' @return `Date` 或 YYYYMMDD 数值向量。
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

#' @rdname num2date
#' @export
date2num <- function(date) {
  as.character(date) %>%
    gsub("-", "", .) %>%
    as.numeric()
}

# 把 YYYYMMDD 数值转为 Date 列
dt_num2date <- function(data) {
  as.data.table(data) %>% mutate(date = num2date(date))
}

# 由 year/month/day 列生成 YYYYMMDD 数值。
ymd2num <- function(data) {
  data$year * 10000 + data$month * 100 + data$day
}

# 由断点位置生成各段的起止索引（含两端）。
segment_bounds <- function(breaks) {
  list(start = c(1, breaks[-length(breaks)] + 1), end = breaks)
}

# 把压缩后的非缺测值放回原位置，其余位置填 NA。
expand_values <- function(x, rows, n) {
  result <- rep(NA_real_, n)
  result[rows] <- round(x, 4)
  result
}

fprintf <- function(fmt, ..., file = "", append = FALSE) {
  cat(sprintf(fmt, ...), file = file, append = append)
}

is.Date <- function(x) is(x, "Date")

plot2 <- function(at, labels, ...) {
  plot(..., xaxt = "n")
  axis(side = 1, at = at, labels = labels)
}
