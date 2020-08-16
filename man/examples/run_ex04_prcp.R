## example 01 ------------------------------------------------------------------
# dat_prcp = read.table("inst/extdata/RHtests_dlyPrcp_ExampleData.txt") %>%
#   data.table() %>% set_names(c("year", "month", "day", "data"))
# use_data(dat_prcp, overwrite = TRUE)

data <- dat_prcp
metadata <- data.table(date = c("1966-11-01", "1976-07-01", "1980-03-01"))
prefix <- "../../OUTPUT/example02/example01"
r <- RHtests_process(data, NULL, metadata, prefix)
plot_RHtests(r)

r <- ReadDLY(dat_prcp)

prefix = "OUTPUT/example04"
U <- FindU.dlyPrcp(dat_prcp, prefix)


# profvis::profvis({
system.time({
  UD <- FindUD.dlyPrcp(dat_prcp, U$TP, prefix)
})
# })
