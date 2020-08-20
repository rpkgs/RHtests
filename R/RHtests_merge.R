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

