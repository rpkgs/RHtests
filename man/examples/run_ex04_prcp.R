## example 01 ------------------------------------------------------------------
# dat_prcp = read.table("inst/extdata/RHtests_dlyPrcp_ExampleData.txt") %>%
#   data.table() %>% set_names(c("year", "month", "day", "data"))
# dat_prcp[data < 0, data:= NA_real_]
# use_data(dat_prcp, overwrite = TRUE)

data <- dat_prcp

metadata <- data.table(date = c("19350927", "19451113", "19460109", "19490113") %>% as.Date("%Y%m%d"))
# prefix <- "../../OUTPUT/example02/example01"
# r <- RHtests_process(data, NULL, metadata, prefix)
# plot_RHtests(r)
r <- ReadDLY(dat_prcp)

prefix = "OUTPUT/example04"
U <- FindU.dlyPrcp(data, prefix)
# if (is_empty(U$TP)) return(NULL)
UD <- FindUD.dlyPrcp(data, U$TP, prefix)
# if (is_empty(UD$TP)) return(NULL)

TP  <- UD$TP
TP2 <- adjust_TP(TP, metadata, maxgap = 366)

# StepSize.dlyPrcp(data, TP2, prefix)
r <- RHtests_stepsize_prcp(data, TP2, prefix)

# system.time({
  # l <- RHtests_process_prcp(data, metadata, prefix)
# })

plot_RHtests(r)

# profvis::profvis({
# system.time({
# })
# })
