test_that("quantile-mapping diagnostics preserve segment-specific adjustments", {
  file <- tempfile()
  on.exit(unlink(file))
  segments <- RHtests:::stepsize_bounds(c(2, 5, 8))
  qm <- list(Nseg.mn = 2, Mq = 2, osmean = matrix(c(1, 2), nrow = 1))
  RHtests:::stepsize_qm_log(qm, segments, 2, 0, function(i) i / 10, file)
  output <- readLines(file)
  expect_true(any(grepl("Adjust to segment 2 : from 3 to 5", output, fixed = TRUE)))
  expect_true(any(grepl("Fcat.Seg1 DP.Seg1 Fcat.Seg2 DP.Seg2", output, fixed = TRUE)))
  expect_true(any(grepl("Seg. 1: mean of QM-adjustments =0.1000", output, fixed = TRUE)))
  expect_false(any(grepl("Seg. 2: mean of QM-adjustments", output, fixed = TRUE)))
  expect_true(any(grepl("Seg. 3: mean of QM-adjustments =0.3000", output, fixed = TRUE)))

  qm$Mq <- 1
  RHtests:::stepsize_qm_log(qm, segments, 2, 0, function(i) stop("not needed"), file)
})
