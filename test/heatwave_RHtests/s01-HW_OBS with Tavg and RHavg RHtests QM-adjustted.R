library(lubridate)
library(RHtests)
library(missInfo)
library(latticeGrob)

load_all()
## 2. 计算热浪特征指标 ---------------------------------------------------------
{
    # file_HI <- "../ChinaHW_data/OUTPUT/OBS_China_HI_met2481-(1951-2018).RDS"
    file_mete_monthly <- "OUTPUT/INPUT_HImete_st2138.rda"
    load(file_mete_monthly)
    df = readRDS("OUTPUT/INPUT_HImete_daily_st2138.RDS")
}

s1_reprocess = FALSE
if (s1_reprocess) {
    ## fill data of date uncontinues values
    df2 <- fix_uncontinue(df)
    # saveRDS(df2, "OUTPUT/INPUT_HImete_daily_st2138.RDS")
}

st = st_met2370[, .(id, site, lon, lat, kind = is_city)]
sites_rural <- st_met2370[is_city == "Rural", site] %>% set_names(., .)
sites = df[, .N, .(site)]$site

# I_bad = c(15, 30, 166, 172, 190)
# I_bad = c(274, 327, 492, 502, 744, 796)
run_RHtests = FALSE
if (run_RHtests){
    InitCluster(10, kill = TRUE)
    # load_all("../meteCMA/RHtests/")
    varname <- "RHavg"
    res  <- RHtests_main(df, sites, varname)
    res2 <- RHtests_rm_empty(res)
    saveRDS(res2, glue("OUTPUT/RHtests_{varname}_monthly_2138.RDS"))

    ## merge yearly and monthly TP
    info  <- tidy_TP(res2)
    info2 <- info[abs(year(date) - year(date_year)) <= 1, ][Idc != "No  ", ]
    sites_adj = info2[, .N, .(site)][, site]
    res_adj = res[sites_adj]
    lst_TP <- split(info2, info2$site)

    out <- RHtests_adj_daily(df, lst_TP, varname)
    saveRDS(out, glue("OUTPUT/RHtests_{varname}_QMadjusted.RDS"))
}

## DEBUG module ----------------------------------------------------------------
{    
    d_Tavg = merge_adjusted(df, "Tavg")
    d_RHavg = merge_adjusted(df, "RHavg")
    df2 = merge(d_Tavg, d_RHavg, by = c("site", "date"))
    setkeyv(df2, c("site", "date"))

    library(JuliaCall)
    julia <- julia_setup()
    julia_source("inst/julia/heat_index.jl")
    df2$Favg = celsius.to.fahrenheit(df2$Tavg) #%>% as.matrix()
    # RH   = as.matrix(df_input$RHavg)
    I_na <- df2[, which(is.na(Tavg + RHavg))]
    # can hold NA values
    HI   <- with(df2[-I_na, ], julia_call("heat_index", Favg, RHavg)) %>% fahrenheit.to.celsius()
    df2$HI <- NA_real_
    df2$HI[-I_na] <- HI
    saveRDS(df2, "OUTPUT/INPUT_HImete_daily_st2138 (QM-adjusted).RDS")
}

## MERGED QM_ADJUST SERIES -----------------------------------------------------
{
    info = df2[, .(n_miss = sum(is.na(HI))), .(site, year = year(date))]
    info$site %<>% as.numeric()
    df_mat <- dcast(df2, date~site, value.var = "HI") 
    r <- HW_index(df_mat, probs = probs)
    r[duration == 0, `:=`(intensity = 0, volume = 0)]
    r2 = info[n_miss < 30, 1:2] %>% merge(r, by = c("site", "year"))
    saveRDS(r2, "../ChinaHW_data/OUTPUT/OBS_China_HW_characteristics_met2481-(1951-2018)-(QM-adjusted).RDS")
}

## DEBUG module ----------------------------------------------------------------
# info = df_month[, .N, .(site, year)][N == 12]
# df_month2 <- merge(df_month, info[, 1:2])
# sitename  <- 53682
{
    # load_all("../meteCMA/RHtests/")
    # sitename    = 53588
    sitename = sites[I_bad[1]]
    d <- df[site == sitename, .(date, Tavg)]
    metadata = st_moveInfo[site == sitename, ] %>% 
        .[period_date_begin > date_begin & 
              period_date_end < date_end, ]
    metadata[, date := period_date_begin]
    
    if (nrow(d) == 0) { message("no data!"); return() }
    l <- RHtests_input(d)
    ## 以monthly为准
    prefix  = "./OUTPUT/example01"
    r_month <- RHtests_process(l$month, NULL, metadata, prefix, is_plot = FALSE, maxgap = 366)
    # r_year  <- RHtests_process(l$year, NULL, metadata, prefix, is_plot = FALSE, maxgap = 366)
}

# res_yearly  = RHtests_main(df_year2 , sites_rural)
# d = df_year2[site == sitename, .(date, Tavg)] %>% format_RHinput()
# r = process_RHtests(d, NULL, metadata)
# res2 = res %>% set_names(sites_rural[1:140]) %>% rm_empty()

# plot_RHtests_multi(res_monthly, "RHtests_monthly.pdf")
# plot_RHtests_multi(res_yearly , "RHtests_yearly.pdf")

## yearly的数据作为辅助验证

# process_RHtests(d, NULL, metadata)
# st_moveInfo[site == 53682, ]
# site moveTimes tag   lon  lat   alt period_date_begin period_date_end date_begin   date_end n_all n_period
# 1: 53682         3   1 11436 3854  94.7        1960-01-01      1964-03-31 1960-01-01 2018-12-31 21550     1552
# 2: 53682         3   2 11441 3838 104.8        1964-04-01      1967-04-30 1960-01-01 2018-12-31 21550     1125
# 3: 53682         3   3 11441 3838 104.1        1967-05-01      2018-12-31 1960-01-01 2018-12-31 21550    18873
# dist
# 1:  0.00000
# 2: 15.46817
# 3: 15.46817

# 现53588站自1998年1月开始站址南迁，海拔高度降低近700m，之后气温明显升高。
# 采用RHtests进行处理非常必要
# st_moveInfo[site == 56584, ]
# site moveTimes tag   lon  lat  alt period_date_begin period_date_end date_begin   date_end n_all n_period
# 1: 56584         3   1 10326 2751 1452        1959-01-01      1978-12-31 1959-01-01 2018-12-31 21915     7305
# 2: 56584         3   2 10315 2742 1452        1979-01-01      2006-12-31 1959-01-01 2018-12-31 21915    10227
# 3: 56584         3   3 10315 2742 1093        2007-01-01      2018-12-31 1959-01-01 2018-12-31 21915     4383
# dist
# 1:  0.00000
# 2: 20.76993
# 3: 20.76993
