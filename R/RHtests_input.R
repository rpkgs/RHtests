#' @importFrom lubridate day
#' @import data.table
#' @export
RHtests_input <- function(d, varnames = NULL) {
  vars_commom <- c("year", "month", "day")
  if (is.null(varnames))
    varnames <- setdiff(colnames(d), c(vars_commom, "date"))

  if (all.equal(intersect(vars_commom, colnames(d)), vars_commom) != TRUE){
      d <- mutate(d, year = year(date), month = month(date), day = day(date))
      # d[, `:=`(year = year(date), month = month(date), day = day(date))]
  }
  d <- d[, .SD, .SDcols = c(vars_commom, varnames)]

  d_year = d[, lapply(.SD, mean, na.rm = TRUE), .SDcols = varnames, .(year)] %>% 
    cbind(month = 1, day = 1) %>%
    reorder_name(vars_commom)

  d_month  = d[, lapply(.SD, mean, na.rm = TRUE), .SDcols = varnames, .(year, month)] %>% 
    cbind(day = 1) %>%
    reorder_name(vars_commom)
  
  list(day = d, month = d_month, year = d_year) %>% lapply(as.data.frame)
}
