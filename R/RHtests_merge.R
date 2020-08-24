merge_TP <- function(res) {
    res2 <- RHtests_rm_empty(res)
    ## merge yearly and monthly TP
    info  <- tidy_TP(res2)
    info2 <- info[abs(year(date) - year(date_year)) <= 1, ][Idc != "No  ", ]
    sites_adj = info2[, .N, .(site)][, site]
    # res_adj = res2[sites_adj]
    lst_TP <- split(info2, info2$site)
    # lst_TP
}

# res2 <- RHtests_doubleCheck(res)
#' rm empty TPs
#' 
#' @param fun intersect or union
#' @export
RHtests_rm_empty <- function(res, fun = intersect) {
    I_left1 <- map(res, 1) %>% which.notnull()
    I_left2 <- map(res, 2) %>% which.notnull()
    I_left <- fun(I_left1, I_left2) # %>% sort()
    names(res[I_left])
    # res[I_left]
}

merge_adjusted <- function(df, varname) {
    infile = glue("OUTPUT/RHtests_{varname}_QMadjusted.RDS")
    out <- readRDS(infile)
    df_adj <- map(out, ~.$data[, .(date, base, QM_adjusted)]) %>% melt_list("site")
    sites_adj <- sort(unique(df_adj$site))

    varnames = c("site", "date", varname)
    df_good = df[!(site %in% sites_adj), .SD, .SDcols = varnames] %>% cbind(QC = 1)
    df_adj2 = df_adj[, .(site, date, QM_adjusted)] %>% set_colnames(varnames) %>% cbind(QC = 0)
    ans = rbind(df_good, df_adj2) %>% set_colnames(c(varnames, paste0("QC_", varname)))
    ans
}


#' @param res2 object returned by [RHtests()]
#' 
tidy_TP <- function(res2) {
    names <- names(res2)
    if (is.null(names)) names <- seq_along(res2)
    sites_TP <- set_names(names, names)

    lst <- foreach(sitename = sites_TP, x = res2, i = icount()) %do% {
        runningId(i, 100)
        # meta = st_moveInfo[site == sitename, ] %>% plyr::mutate(date = period_date_begin)
        # meta$date %<>% as.Date()

        year = x$year$TP
        month = x$month$TP
        # year$date %<>% as.character() %>% as.Date("%Y%m%d")
        # month$date %<>% as.character() %>% as.Date("%Y%m%d")

        # for(j in 1:nrow(month)) {
        month2 <- foreach(j = 1:nrow(month)) %do% {
            date = month$date[j]
            diff_year = difftime(date, year$date, units = "days") %>% as.numeric()
            I_year = which.min(abs(diff_year))

            diff_meta = difftime(date, month$date_meta, units = "days") %>% as.numeric()
            I_meta = which.min(abs(diff_meta))

            cbind(month[j, 1:9],
                date_year = year$date[I_year], day2_year = diff_year[I_year],
                date_meta = month$date_meta[I_meta], day2_meta = diff_meta[I_meta]
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
