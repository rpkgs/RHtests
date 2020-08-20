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

#' Revise site moving information for all sites
#'
#' @param info site moving info returned by [get_moveInfo()]
#' @param dist_max If the site moving distance beyond the `dist_max`, it will
#' be regarded as outlier and will be fixed by [revise_location()].
#'
#' @details
#' Continuous moving location outliers are also considered at here.
#'
#' @export
revise_locaiton_main <- function(info, dist_max = 50) {
    info$QC <- ""
    sites <- info[dist >= dist_max]$site %>% unique()
    info_bad <- info[site %in% sites, ]

    temp <- foreach(sitename = sites, i = icount()) %do% {
        d <- info[site == sitename, ]
        # d$QC = "raw"
        prefix <- sprintf("%02dth:%s", i, sitename)
        r <- revise_locaiton(d, prefix, dist_max)
        # r$status
        r
    }
    # %>% reorder_name(c("site", "moveTimes", "tag", "lon", "lat", "alt", "dist", "n_period", "QC"))
    df <- do.call(rbind, temp)

    d_fixed <- df[QC %in% c("good", "margin"), .N, .(site)] # %>% nrow()
    d_unsure <- df[QC %in% c("suspicious"), .N, .(site)] # %>% nrow()
    d_unfixed <- df[QC %in% c("unfixed"), .N, .(site)]
    sites_bad <- d_unsure$site
    sites_unfixed <- setdiff(d_unfixed$site, d_unsure$site)
    # n_fixed = length(sites) - length(sites_bad)
    ok(sprintf(
        "[info] %d sites fixed, %d sites unfixed, %s sites not sure\n",
        nrow(d_fixed), length(sites_unfixed), length(sites_bad)
    ))

    info_final <- rbind(info[!(site %in% sites)], df)
    info_final
    # sites_bad = sites[which(unlist(temp) == "bad")]
    # sites_bad = c("51058", "52378", "52607", "52884", "53730", "54287")
}
