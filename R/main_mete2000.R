#' fix un-continuous date of the INPUT data
#' @param df A data.frame with the columns of `site`, `date` and others
#' @export
fix_uncontinue <- function(df) {
    df2 = ddply(df, .(site), function(d) {
        n = nrow(d)
        date_begin = d$date[1]
        date_end   = d$date[n]
        if (month(date_begin) != 1) 
            date_begin = make_date(year(date_begin) + 1, 1, 1)
        if (month(date_end) != 12) 
            date_end = make_date(year(date_end) - 1, 12, 31)
        
        temp = data.table(site = d$site[1], 
                          date = seq.Date(date_begin, date_end, by = "day"))
        merge(d, temp, c("site", "date"), all.y = TRUE)
    }, .progress = "text")
    df2
}
