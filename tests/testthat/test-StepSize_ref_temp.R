test_that("RHtests_process temp works", {
  # expect_equal(2 * 2, 4)
  l <- RHtests_input(dat_temp) #%>% str()
  # prefix <- "../../OUTPUT/example02/example02"
  prefix <- "OUTPUT/example02"
  B <- l$month[, c(1, 2, 3, 4)]
  R <- l$month[, c(1, 2, 3, 5)]

  r <- RHtests_process(B, R, metadata = NULL, prefix)
  plot_output(r$data)

  expect_equal(nrow(r$data), 744)
  expect_equal(r$TP$date %>% date2num(), c(198812, 200606, 200910) * 100 + 1)
})

# px <- proffer::pprof({
#   r <- RHtests_process(B, R, metadata = NULL, prefix, is_plot = FALSE)
# })
