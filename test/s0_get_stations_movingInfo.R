library(plyr)
library(glue)
library(foreach)
library(iterators)
library(lubridate)
library(missInfo)
library(purrr)

# dirs = list.dirs("N:/DATA/China/2400climate data")
dir_root = "N:/DATA/China/2400climate data" %>% path.mnt()
varnames = c("EVP", "GST", "PRE", "PRS", "RHU", "SSD", "TEM", "WIN")

# [826] 69473, 59034 2429 10839   1851 2019 10 235000002     64      5    106      5 9 9 9 9 9
#   win_avg: 5000002 -> 20
# [827] 17129: 53955 3531 11028   4616 2019 11 292500012     33      3     64      3 9 9 9 9 9
#   win_avg: 500012 -> 12
# [829]  3792: 51238 4454  8204   5322 2020  1  72500001     21     13     32     13 9 9 9 9 9
#   win_avg: 2500001 -> 10
foreach(varname = varnames[3], i = icount()) %do% {
    outfile = glue("{dir_root}/SURF_CLI_CHN_MUL_DAY_{varname} (195101-202003).csv")
    # if (file.exists(outfile)) return()

    indir = glue("{dir_root}/{varname}")
    files = dir(indir, "*.TXT", full.names = TRUE)

    lst <- foreach(file = files[1:length(files)], i = icount()) %do% {
        runningId(i, 10)
        tryCatch({
            fread(file)
        }, warning = function(e) {
            message(sprintf('[i] %s: %s', i, basename(file), e$message))
        })
    }
    # lst <- llply(files, fread, .progress = "text")
    df = do.call(rbind, lst)
    invisible()
    # fwrite(df, outfile)
}

# ------------------------------------------------------------------------------
dir_root <- "N:/DATA/China/2400climate data" %>% path.mnt()
varnames <- c("EVP", "GST", "PRE", "PRS", "RHU", "SSD", "TEM", "WIN")
vars_common = c("site", "lat", "lon", "alt", "year", "month", "day")

files = dir(dir_root, "*.csv", full.names = TRUE) %>% set_names(varnames)
# I_sel = 1:7
lst <- map(files[7], ~fread(.x, select = I_sel) %>% set_colnames(vars_common[I_sel]))
l   <- map(lst, ~.x[year == 2019 & month == 12]) %>% map_int(nrow)

saveRDS(lst, "a.RDS", compress = FALSE)
# 2019年数据异常
# ds$WIN[N < 30]
# site  N
# 1: 51468 11
# 2: 54287 20
# 3: 54646 20
# 4: 57713 19
# 5: 58358 20
# 6: 58726 22

# 453
# d = fread(files[1])
## 获取台站变迁记录
# get_history_location <- function() {
df = lst$TEM
st_full      = df[, .(site, lat, lon, alt, date = make_date(year, month, day))]
st_1961_2018 = st[date >= "1960-01-01" & date <= "2018-12-31"]
# st = st[date <= "2019-12-31"]
info_full = get_moveInfo(st_full)

{
    load_all("../missInfo/")
    # 1961_2018, 132 unchanged: lon^2 + lat^2 + alt
    info_1961_2018 = get_moveInfo(st_1961_2018)
    # 337 unchanged: lon^2 + lat^2

}

# 1. 考虑连续型的bad values

{
    dist_max = 50
    info_1961_2018$QC = ""
    sites = info_1961_2018[dist >= 50]$site %>% unique()
    info = info_1961_2018[site %in% sites, ]

    info = info_1961_2018[site %in% sites, ]
    temp = foreach(sitename = sites, i = icount()) %do% {
        d = info[site == sitename, ]
        # d$QC = "raw"
        prefix = sprintf("%02dth:%s", i, sitename)
        r <- revise_locaiton(d, prefix, dist_max)
        # r$status
        r
        # print(d)
    }
    # %>% reorder_name(c("site", "moveTimes", "tag", "lon", "lat", "alt", "dist", "n_period", "QC"))
    df <- do.call(rbind, temp)

    d_fixed   = df[QC %in% c("good", "margin"), .N, .(site)] #%>% nrow()
    d_unsure  = df[QC %in% c("suspicious"), .N, .(site)] #%>% nrow()
    d_unfixed = df[QC %in% c("unfixed"), .N, .(site)]
    sites_bad = d_unsure$site
    sites_unfixed = setdiff(d_unfixed$site, d_unsure$site)
    # n_fixed = length(sites) - length(sites_bad)
    ok(sprintf("[info] %d sites fixed, %d sites unfixed, %s sites not sure\n",
               nrow(d_fixed), length(sites_unfixed), length(sites_bad)))
    info_final = rbind(info_1961_2018[!(site %in% sites)], df)
    # sites_bad = sites[which(unlist(temp) == "bad")]
    # sites_bad = c("51058", "52378", "52607", "52884", "53730", "54287")
}

## 修复了26/51个站点

# use_data(st_moveInfo, overwrite = TRUE)
# 89站点
# }
## add moving dist at here

