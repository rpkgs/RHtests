## example 02 ------------------------------------------------------------------
Bseries <- system.file("extdata/Example2.dat", package = "RHtests")
Rseries <- system.file("extdata/Example2_Ref.dat", package = "RHtests")

metadata <- data.table(date = c("19740200", "19751100"))

r <- RHtests_process(Bseries, Rseries, metadata, outdir = "OUTPUT/example02")
plot_RHtests(r, "RHtests-ref.pdf")

# expect_equal(nrow(r$turningPoint), 11)
# expect_equal(as.character(r$turningPoint[kind == 0]$date),
#               metadata$date)
# expect_equal(nrow(r$data), 612)
