test_that("point tables respect append and preserve their headers", {
  point <- data.table::data.table(
    kind = 1, Idc = "Yes ", date = 20000101, Ic = 12, Nseg = 24,
    stepsize = 1.23456, probL = 0.8, probU = 0.95, plev = 0.95,
    PFx = 5, PFx95l = 3, PFx95h = 6, prob = 0.9, PFx95 = 4
  )
  file <- tempfile()
  on.exit(unlink(file))

  for (reference in c(FALSE, TRUE)) {
    writeLines("old content", file)
    result <- RHtests:::write_points_table(list(point), file,
      reference = reference, stepsizes = 1.23456
    )
    lines <- readLines(file)
    expect_false(any(grepl("old content", lines, fixed = TRUE)))
    if (!reference) expect_match(lines[1], "^#  1 changepoints")
    header <- if (reference) 1 else 2
    expect_match(lines[header], "^kind")
    expect_equal(result$stepsize, 1.2346)

    writeLines("previous content", file)
    RHtests:::write_points_table(list(point), file,
      reference = reference, append = TRUE, stepsizes = 1.23456
    )
    lines <- readLines(file)
    expect_identical(lines[1], "previous content")
    expect_match(lines[header + 1], "^kind")
  }

  writeLines("old content", file)
  expect_length(RHtests:::write_points_table(list(), file), 0)
  expect_length(readLines(file), 0)
  writeLines("previous content", file)
  expect_length(RHtests:::write_points_table(list(), file, append = TRUE), 0)
  expect_identical(readLines(file), "previous content")
})
