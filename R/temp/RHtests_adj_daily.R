#' get the daily QM-adjusted results
#' @export
RHtests_adj_daily <- function(df, lst_TP, varname = "Tavg") {
    sites <- names(lst_TP) %>% set_names(., .)
    res_daily <- foreach(sitename = sites, TP = lst_TP, i = icount()) %dopar% {
        runningId(i)
        tryCatch({
            # d <- df[site == sitename, .(date, Tavg)] %>% format_RHinput()
            d <- df[site == sitename, .SD, .SDcols = c("date", varname)]
            metadata = st_moveInfo[site == sitename, ]
            metadata[, date := period_date_begin]

            if (nrow(d) == 0) { message("no data!"); return() }
            l <- RHtests_input(d)
            prefix  = "./OUTPUT/example01"
            # r_month <- RHtests_process(l$month, NULL, metadata, prefix, is_plot = FALSE, maxgap = 366)
            # r_year  <- RHtests_process(l$year, NULL, metadata, prefix, is_plot = FALSE, maxgap = 366)
            ## need to input TP
            # TP <- r_month$TP
            r_daily <- RHtests_stepsize(l$day, NULL, TP, prefix = prefix, is_plot = FALSE)
            # list(year = r_year, month = r_month)
        }, error = function(e) {
            message(sprintf("%s", e$message))
        })
    }
    res_daily
}
