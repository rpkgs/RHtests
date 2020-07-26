## example 01 ------------------------------------------------------------------
data <- dat_example01
metadata <- data.table(date = c("1966-11-01", "1976-07-01", "1980-03-01"))
prefix <- "../../OUTPUT/example02/example01"
r <- RHtests_process(data, NULL, metadata, prefix)
plot_RHtests(r)
