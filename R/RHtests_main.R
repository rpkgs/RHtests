# num2date <- function(x) {
#     as.character(x) %>% gsub("00$", "01", .) %>% as.Date("%Y%m%d")
# }

format_RHinput <- function(d) {
    varnames <- setdiff(colnames(d), c("site", "date"))
    d[, .(year = year(date), month = month(date), day = day(date))] %>%
        cbind(d[, ..varnames])
}

#' get the turning points of monthly and yearly data
#'
#' @param df A data.frame with the columns at least of `site`, 'date', 'varname'
#' @param st_moveInfo
#'
#' @export
RHtests_main <- function(df, st_moveInfo, sites, varname, .parallel = FALSE) {
    sites %<>% set_names(., .)
    res <- foreach(sitename = sites, i = icount()) %dopar% {
        # for(i in seq_along(sites_rural[1:30])) {
        runningId(i)
        # sitename = sites_rural[i]
        tryCatch({
            # d <- df[site == sitename, .(date, Tavg)] %>% format_RHinput()
            d <- df[site == sitename, .SD, .SDcols = c("date", varname)]
            date_begin = d$date[1]
            date_end = d$date[nrow(d)]
            metadata = st_moveInfo[site == sitename, ] %>%
                .[period_date_begin > date_begin &
                    period_date_end < date_end, ]
            metadata[, date := period_date_begin]

            if (nrow(d) == 0) { message("no data!"); return() }
            l <- RHtests_input(d)
            ## 以monthly为准
            prefix  = "./OUTPUT/example01"
            # browser()
            r_month <- RHtests_process(l$month, NULL, metadata, prefix, is_plot = FALSE, maxgap = 366)
            r_year  <- RHtests_process(l$year, NULL, metadata, prefix, is_plot = FALSE, maxgap = 366)
            # TP <- r_month$TP
            # r_daily <- RHtests_stepsize(l$day, NULL, TP, prefix = prefix, is_plot = TRUE)
            list(year = r_year, month = r_month)
        }, error = function(e) {
            message(sprintf("[%d] %s", i, e$message))
        })
    }
    res
    # temp <- map(res, 2) %>% rm_empty()
    # sites_nonhomo <- names(temp)
    # TP <- melt_list(temp, "site")
    # res2 <- res[sites_nonhomo]
    # listk(sites_nonhomo, TP, result = res2)
}

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


tidy_TP <- function(res2) {
    sites_TP <- names(res2) %>% set_names(., .)
    lst <- foreach(sitename = sites_TP, x = res2, i = icount()) %do% {
        runningId(i, 100)
        meta = st_moveInfo[site == sitename, ] %>% plyr::mutate(date = period_date_begin)
        meta$date %<>% as.Date()

        year = x$year$TP
        month = x$month$TP
        # year$date %<>% as.character() %>% as.Date("%Y%m%d")
        # month$date %<>% as.character() %>% as.Date("%Y%m%d")

        # for(j in 1:nrow(month)) {
        month2 <- foreach(j = 1:nrow(month)) %do% {
            date = month$date[j]
            diff_year = difftime(date, year$date, units = "days") %>% as.numeric()
            I_year = which.min(abs(diff_year))

            diff_meta = difftime(date, meta$date, units = "days") %>% as.numeric()
            I_meta = which.min(abs(diff_meta))

            cbind(month[j, 1:9],
                date_year = year$date[I_year], day2_year = diff_year[I_year],
                date_meta = meta$date[I_meta], day2_meta = diff_meta[I_meta]
            ) %>%
                reorder_name(c(
                    "kind", "Idc", "Ic", "date", "date_meta", "date_year",
                    "day2_meta", "day2_year"))
        } %>% do.call(rbind, .)
        cbind(site = sitename, month2)
        # listk(TP = month2, meta)
    }
    info <- do.call(rbind, lst)
    info
}

plot_RHtests_multi <- function(obj, outfile = "RHtests.pdf") {
    dout <- map(obj$result, ~ .$data[, .(date = num2date(date), base, QM_adjusted)]) %>%
        melt_list("site")
    n = length(obj$result)

    p <- ggplot(dout, aes(date, y = QM_adjusted - base)) +
        geom_line() +
        # geom_line(aes(date, QM_adjusted), color = "blue") +
        facet_wrap(~site, scales = "free", ncol = 2)
    write_fig(p, outfile, 10, 50/70*n)
}
