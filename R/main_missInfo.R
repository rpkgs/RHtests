#' The distance of locations from the location having the
#'
#' @export
get_dist <- function(lon, lat, n_period) {
    P = cbind(lon, lat) %>% deg2dec()
    i = which.max(n_period) #

    if (nrow(P) > 1) {
        dist = rdist.earth(P[i,,drop = FALSE], P[,,drop = FALSE])[1, ]
        round(dist, 2)
    } else 0
}

#' guess the exact location based on distance score
dist_score <- function(x, y) {
    (x$lon == y$lon) + (x$lat == y$lat) + (abs(x$alt - y$alt) <= 0.005)
}

last <- function(x) {x[length(x)]}

#' Revise meteorological station's location error
#'
#' @param d A data.frame with the columns at least of `dist`, `lon`, `lat`, `alt`
#' and `n_period`, `QC`
#' @param vars_rm variables not shown in the console
#' @param score_min if distance score lower than `score_min`, will not fix
#'
#' @return
#' sites with not outliers, QC = "";
#' 
#' sites with not outliers, for outliers, QC equals:
#' score = 0: unfixed
#' score = 1: marginal quality fixed
#' score = 2: good quality fixed
#' 
#' @export
revise_locaiton <- function(d, prefix = "",
    dist_max = 100,
    score_min = 2,
    vars_rm = c("date_begin", "date_end", "n_all"),
    verbose = 2)
{
    vars = setdiff(colnames(d), vars_rm)
    n = length(d$dist)
    # grps = cumsum(1, diff(d$dist >= dist_max) != 0)
    inds = which(d$dist >= dist_max)
    grps = if (length(inds) == 1) {
        cumsum(c(1, diff(inds) != 1))
    } else 1
    ngrp = max(grps)

    num = 0
    for (i in 1:ngrp) {
        ind = inds[grps == i]

        i_prev = ind[1] - 1
        if (i_prev < 1) i_prev <- NULL
        i_next = ind[length(ind)] + 1
        if (i_next > n) i_next <- NULL
        ind_inspect = c(i_prev, ind, i_next)

        i_candinate = c(i_prev, i_next)

        # print(i_candinate)
        for (j in ind) {
            dist = d$dist[j]
            scores = dist_score(d[i_candinate, ], d[j, ])
            i_opt = i_candinate[which.max(scores) %>% last()]
            score = max(scores)
            # browser()
            # if (dist >= 500), we have to fix it
            is_good = score >= 2 || (dist > 500 && score >= 1) # confidential fix
            is_print = verbose == 1 || (verbose >= 2 && !is_good)
            if (is_print) {
                num = num + 1
                if (num == 1)
                    ok(sprintf("[%s] ==============================================================\n", prefix))
                fprintf("---------------------------------------------------------\n")
                print(d[ind_inspect, ..vars])
                fprintf("---\n")
            }

            is_low_score = FALSE
            if (is_good) {
                d$lat[j] = d$lat[i_opt]
                d$lon[j] = d$lon[i_opt]
                d$alt[j] = d$alt[i_opt]
                d$QC[j]  = "good"
            } else {
                warn(sprintf("[w]: dist_score < %d\n", score_min))
                if (score >= score_min) {
                    d$lat[j] = d$lat[i_opt]
                    d$lon[j] = d$lon[i_opt]
                    d$alt[j] = d$alt[i_opt]
                    d$QC[j]  = "marginal"
                } else {
                    d$QC[j] = ifelse(score != 0, "suspicious", "unfixed")
                    warn(sprintf("[w]: not fixed, too low score=%d\n", min(scores)))
                    is_low_score = TRUE
                }
            }
            d[, dist := get_dist(lon, lat, n_period)]
            if (is_print && !is_low_score) print(d[ind_inspect, ..vars])
        }
    }
    d
    # listk(d, status)
}

get_moveInfo <- function(st, prefix = NULL) {
    st$lon %<>% deg2dec() %>% dec2deg()
    st$lat %<>% deg2dec() %>% dec2deg()
    st$alt %<>% get_alt()

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
        fwrite(st_moveInfo, glue::glue("data-raw/mete2481_站点变迁记录-{prefix}-({str_begin}-{str_end}).csv"))
    }
    st_moveInfo
}

