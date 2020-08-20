
st_moveInfo <- function(st, prefix = NULL) {
    st$lon %<>% deg2dec() %>% dec2deg()
    st$lat %<>% deg2dec() %>% dec2deg()
    st$alt %<>% get_alt()
    n = length(unique(st$site)) # length of sites

    st_moveInfo = ddply(st, .(site), function(d) {
        # d = st[site == 58246]
         # + alt
        d$tag = d[, lon^2 + lat^2] %>%
            {c(1, abs(diff(.)) > 0.1)} %>% cumsum()
        # d$date = d[, make_date(year, month, day)]
        date_begin = min(d$date)
        date_end   = max(d$date)

        d[, .(period_date_begin = min(date), period_date_end = max(date),
              date_begin, date_end),
            .(site, tag, lon, lat, alt)]
    }, .progress = "text")

    # st_moveInfo %<>% do.call(rbind, .)
    st_moveInfo[, moveTimes := max(tag), .(site)]
    st_moveInfo %<>% reorder_name(c("site", "moveTimes", "tag"))
    st_moveInfo[, 7:10] %<>% map(as.Date)
    # st_moveInfo[, alt := get_alt(alt)]
    st_moveInfo[, `:=`(n_all = difftime(date_end, date_begin) %>% as.numeric() %>% add(1),
                       n_period = difftime(period_date_end, period_date_begin, units = "days") %>% as.numeric() %>% add(1))]

    st_moveInfo[, dist := get_dist(lon, lat, n_period), .(site)]
    if (!is.null(prefix)) {
        str_begin = st[1, format(date, "%Y%m")]
        str_end   = st[nrow(st), format(date, "%Y%m")]
        fwrite(st_moveInfo, glue::glue("data-raw/mete{n}_站点变迁记录-{prefix}-({str_begin}-{str_end}).csv"))
    }
    st_moveInfo
}

#' @export
retidy_moveInfo <- function(st) {
    st_moveInfo = ddply(st, .(site), function(d) {
        # d = st[site == 58246]
         # + alt
        d$tag = d[, lon^2 + lat^2] %>%
            {c(1, abs(diff(.)) > 0.1)} %>% cumsum()
        # d$date = d[, make_date(year, month, day)]
        date_begin0 = d$date_begin[1]
        date_end0   = d$date_end[1]

        d[, .(period_date_begin = min(period_date_begin), 
              period_date_end = max(period_date_end),
              date_begin = date_begin0, date_end = date_end0),
            .(site, tag, lon, lat, alt)]
    }, .progress = "text")
    st_moveInfo
}
