
# 948,
sites_worst = names(lst_TP)         # monthly and yearly all confirmed
sites_worse = RHtests_rm_empty(res) # monthly and yearly all have TPs
sites_bad   = RHtests_rm_empty(res, union)

sites_good = setdiff(sites, sites_bad) # 410

info =  info2_1961_2018[site %in% sites_good, .(alt_diff = diff(range(alt))), .(site)]

sites_unchanged = info2_1961_2018[, .N, .(site)][N == 1]
# only 130 sites location unchanged
