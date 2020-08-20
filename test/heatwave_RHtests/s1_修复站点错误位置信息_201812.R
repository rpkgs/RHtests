st_moveInfo <- fread("data-raw/mete2481_站点变迁记录-(195101-201812).csv")

# 1961年之前，数据有较多的错误
st = st_moveInfo[period_date_end > "1951-01-01", ]
st[, dist := get_dist(lon, lat), .(site)]
info = st[dist > 100, .N, .(site)][order(-N), ]

# 2019年数据错误较多

# dist in the unit of km
st = st_moveInfo[period_date_end > "1961-01-01", ]
st[, dist := get_dist(lon, lat), .(site)]
info = st[dist > 100, .N, .(site)][order(-N), ]
# sites_bad = info$site
sites_bad_from1951 = c(53938, 54830, 51573, 51814, 52378, 53693, 50915, 57067, 53195, 53730,
                       51058, 56193, 57761, 51329, 51722, 52607, 52884, 54287, 54313, 55279, 55294, 55578, 56792, 57080, 57460, 58246)
sites_bad_from1961 = c(52378, 53730, 51058, 51329, 51722, 52607, 52884, 53195, 54287, 54313,
                       55279, 55294, 55578, 56792, 57080, 57460, 58246)
# st[site == sites_bad_from1961[1]]
all.equal(info$site, sites_bad_from1961)

st_bad = st[site %in% sites_bad_from1961]
# fwrite(st_bad, "data-raw/st_195101-201812-origin.csv")

st_patch <- fread("data-raw/st_195101-201812-patches.csv")
# st_patch$alt %<>% multiply_by(100)
# fwrite(st_patch, "data-raw/st_195101-201812-patches.csv")
st_patch[,  dist := get_dist(lon, lat), .(site)]
st_patch[dist > 100]
st_patch = reflag(st_patch)

st2 <- st[!(site %in% sites_bad_from1961), ] %>% rbind(st_patch)
setkeyv(st2, c("site", "tag"))

st_moveInfo <- st2
use_data(st_moveInfo, overwrite = TRUE)
fwrite(st_moveInfo, "data-raw/st_195101-201812-adjusted.csv")

# sites_bad2 = st2[dist > 50, site] %>% unique()
# st2[site %in% sites_bad2, ]
# sitenames <- get_sitenames(st_met2481$site[1:10])
