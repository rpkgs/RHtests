test_that("StepSize with noref works", {
  data <- dat_example01
  metadata <- data.table(date = c("1966-11-01", "1976-07-01", "1980-03-01"))

  prefix <- "../../OUTPUT/example02/example01"
  r <- RHtests_process(data, NULL, metadata, prefix, maxgap = 366)
  plot_RHtests(r)

  expect_equal(nrow(r$data), 636)
  expect_equal(r$TP$date %>% date2num(), c(197001, 197402, 197607) * 100 + 1)
})
