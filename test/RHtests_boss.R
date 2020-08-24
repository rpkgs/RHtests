library(plyr)
library(glue)
library(foreach)
library(iterators)
library(lubridate)
library(missInfo)
library(purrr)
library(matrixStats)

get_refer_sites <- function(df2, varname = "Tavg") {
    date = seq(ymd("1951-01-01"), ymd("2019-12-31"), by = "day")
    mat = dcast(df2, date ~ site, value.var = varname)[, -1] %>% as.matrix()

    ## when aggregate daily to monthly scale, if more than 3 invalid values, monthly
    # value will be set to NA
    mat_month = apply_col(mat, by = format(date, "%Y-%m-01"))
    mat_month_miss = apply_col(is.na(mat), by = format(date, "%Y-%m-01"), colSums2)
    mat_month[mat_month_miss >= 3] = NA_real_

    ## searching potential reference sites
    # check site names order first
    if (!all.equal(st$site %>% as.character(), colnames(mat)[-1])) {
        stop("check site names order first!")
    }

    # 1. 站点长度至少大于30年（76个站点被移除）
    st_refs = st_refer(st, dist, mat_month)
    st_refs_opt = st_refer_opt(st_refs, sites_worst)
    d_refs <- melt_list(st_refs_opt %>% rm_empty(), "target")
    # save(st_refs, st_refs_opt, d_refs, file = "OUTPUT/mete2481_TEM_st_refer.rda")
}

get_metadata <- function(d, sitename, st_moveInfo) {
    date_begin <- d$date[1]
    date_end <- d$date[nrow(d)]
    metadata <- st_moveInfo[site == sitename, ] %>%
        .[period_date_begin > date_begin & period_date_end < date_end, ] %>%
        .[, date := period_date_begin]
}

# ------------------------------------------------------------------------------
dir_root <- "N:/DATA/China/2400climate data" %>% path.mnt()
varnames <- c("EVP", "GST", "PRE", "PRS", "RHU", "SSD", "TEM", "WIN")
# merge_mete2000_txts(dir_root, is_save = TRUE)
vars_common = c("site", "lat", "lon", "alt", "year", "month", "day")

# 2019-2020年数据异常，观测数据不是完整的一整月的数据，
# df = lst$TEM
file_met = "OUTPUT/mete2481_Tavg_daily (195101-202003).rda"
if (!file.exists(file_met)) {
    files = dir(dir_root, "*.csv", full.names = TRUE) %>% set_names(varnames)
    I_sel = 1:13
    read_data <- function(file) fread(file, select = I_sel) %>% set_colnames(vars_common[I_sel])
    # lst <- map(files[7], read_data)
    df <- read_data(files[7])
    obs_types <- c("avg", "max", "min")
    prefix    <- "T"
    varnames  <- c(paste0(prefix, obs_types), paste0("QC_", prefix, obs_types))
    df[df == 32766] = NA_integer_
    df[, date := make_date(year, month, day)]

    st_full      = df[, .(site, lat, lon, alt, date = make_date(year, month, day))]
    st_1961_2018 = st_full[date >= "1960-01-01" & date <= "2018-12-31"]

    # load_all("../missInfo/")
    # 1961_2018, 132 unchanged: lon^2 + lat^2 + alt
    info_1961_2018  = st_moveInfo(st_1961_2018)
    info_full       = st_moveInfo(st_full)
    info2_1961_2018 = revise_locaiton_main(info_1961_2018, dist_max =  50)
    info2_full      = revise_locaiton_main(info_full, dist_max =  50)
    # 337 unchanged: lon^2 + lat^2

    st = ddply(info2_1961_2018, .(site), . %>% .[which.max(n_period), ])
    coord = st[, .(lon = deg2dec(lon), lat = deg2dec(lat))] %>% as.matrix()
    dist = rdist.earth(coord)

    df2 = df[, .(site, date, Tavg)] %>% fix_uncontinue(complete_year = TRUE)
    save(df2, info2_1961_2018, info2_full, st, dist, file = file_met)
} else {
    load(file_met)
}
info2_1961_2018 <- retidy_moveInfo(info2_1961_2018)

## 修复了26/51个站点
# use_data(st_moveInfo, overwrite = TRUE)
# nrow(df)/nrow(df2)

varname = "Tavg"
sites = st$site
nsite = length(sites)
file_RHtests_monthly = glue("OUTPUT/RHtests_mete{nsite}_{varname}_monthly.RDS")
file_RHtests_daily   = glue("OUTPUT/RHtests_{varname}_QMadjusted.RDS")

InitCluster(14)
## 1. PMF test -----------------------------------------------------------------
res  <- RHtests_main(df2, st_moveInfo = info2_1961_2018, sites, varname)
saveRDS(res, file_RHtests_monthly)

res = readRDS(file_RHtests_monthly)
lst_TP <- merge_TP(res)
out <- RHtests_adj_daily(df2, lst_TP, varname)
saveRDS(out, file_RHtests_daily)

# monthly和yearly一致的，884 TPs
# temp <- RHtests_main(df2, st_moveInfo = info2_1961_2018, sitename, varname)

st_moveInfo = info2_1961_2018
varname = "Tavg"
d_refs = get_refer_sites(df2, varname)
sites_miss = setdiff(sites, d_refs$target) %>% as.character()

# 剩余的153个站点中，64个存在显著的突变点，21个不存在突变点，其余monthly和yearly数据的突变点不一致
#  81 not ref-sites
# 153 not ref-sites, sites_worst
# 308 not ref-sites, sites_worse
lst <- split(df2[, -1], df2$site)
# setkeyv(df2, c("site", "date"))

for (i in 1:nrow(d_refs)) {
    site_target = d_refs$target[i]
    site_refer  = d_refs$site[i]
    i_t = match(site_target, sites)
    i_r = match(site_refer, sites)

    d_target = lst[[i_t]]
    d_refer  = lst[[i_r]]
    d = merge(d_target, d_refer %>% set_names(c("date", "ref")), all.x = TRUE)
    metadata = get_metadata(d, site_target, st_moveInfo)

        
}
