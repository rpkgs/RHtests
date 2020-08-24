test_that("StepSize.wRef works", {
  Bseries <- system.file("extdata/Example2.dat", package = "RHtests")
  Rseries <- system.file("extdata/Example2_Ref.dat", package = "RHtests")

  metadata <- data.table(date = c("1974-02-01", "1975-11-01"))
  prefix <- "../../OUTPUT/example02/example02"

  # r <- RHtests_process(Bseries, Rseries, metadata, prefix)
  # plot_RHtests(r)

  # Test the impact of NA values on RHtests
  B = read.table(Bseries, na.strings = "-999.99") %>% data.table()
  R = read.table(Rseries, na.strings = "-999.99") %>% data.table()
  R[V1 >= 2002, V4:=NA]

  r <- RHtests_process(B, R, metadata, prefix)
  plot_RHtests(r)

  # expect_equal(nrow(r$TP), 11)
  expect_equal(as.character(r$TP[kind == 0]$date),
               metadata$date)
  expect_equal(nrow(r$data), 612)
})
