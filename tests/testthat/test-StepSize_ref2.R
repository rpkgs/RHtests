test_that("multiplication works", {
    # expect_equal(2 * 2, 4)
  l <- RHtests_input(dat_temp) #%>% str()
  prefix <- "../../OUTPUT/example02/example02"

  prefix <- "OUTPUT/example02/example02"
  B <- l$month[, c(1, 2, 3, 4)]
  R <- l$month[, c(1, 2, 3, 5)]

  r <- RHtests_process(B, R, metadata = NULL, prefix)
  plot_RHtests(r)


})
