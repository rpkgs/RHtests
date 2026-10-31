test_that("fprintf writes formatted output to files and standard output", {
  file <- tempfile()
  on.exit(unlink(file))
  RHtests:::fprintf("value: %d\n", 1, file = file)
  RHtests:::fprintf("value: %d\n", 2, file = file, append = TRUE)
  expect_equal(readLines(file), c("value: 1", "value: 2"))
  expect_output(RHtests:::fprintf("value: %d", 3), "value: 3")
})
