test_that("changepoint logs use the corresponding test statistic", {
  result <- data.table::data.table(
    Ic = 12, date = 20000101, kind = 1, prob = 0.9,
    probL = 0.8, probU = 0.95, PFx = 5, PFx95 = 4,
    PFx95l = 3, PFx95h = 6, Nseg = 24
  )
  file <- tempfile()
  on.exit(unlink(file))
  RHtests:::write_test_log(list(result), "PMF", file)
  RHtests:::write_test_log(list(result), "PMT", file)
  lines <- readLines(file)
  expect_match(lines[1], "^PMF :.*PFmax= *5\\.0000")
  expect_match(lines[2], "^PMT :.*PTmax= *5\\.0000")
  expect_match(lines[1], "CV95= *4\\.0000\\( *3\\.0000- *6\\.0000\\)")
})
