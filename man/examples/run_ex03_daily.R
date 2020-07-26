
load("debug.rda")
prefix <- "OUTPUT/example03/example03"

l <- RHtests_input(d)

## 以monthly为准
r_month <- RHtests_process(l$month, NULL, metadata, prefix)
r_year  <- RHtests_process(l$year, NULL, metadata, prefix)

TP <- r_month$TP
r_daily <- RHtests_stepsize(l$day, NULL, TP, prefix = prefix, is_plot = TRUE)

# temp <- Read(d)
# U <- FindU(output = prefix, is_plot = is_plot)
#
# profvis::profvis({
#   U <- FindU(output = prefix, is_plot = is_plot)
#   # r <- process_RHtests(d, NULL, metadata, prefix)
# })
