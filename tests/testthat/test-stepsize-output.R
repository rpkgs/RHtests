test_that("StepSize output and plotting do not change computed results", {
  points <- data.table::data.table(date = 19740201, kind = 1)
  base <- system.file("extdata/Example1.dat", package = "RHtests")
  base_ref <- system.file("extdata/Example2.dat", package = "RHtests")
  ref <- system.file("extdata/Example2_Ref.dat", package = "RHtests")
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  for (reference in c(FALSE, TRUE)) {
    run <- function(prefix, debug, is_plot) {
      if (reference) {
        StepSize.wRef(base_ref, ref, points, prefix, debug = debug, is_plot = is_plot)
      } else {
        StepSize(base, points, prefix, debug = debug, is_plot = is_plot)
      }
    }
    plain <- file.path(dir, paste0("plain", reference))
    report <- file.path(dir, paste0("report", reference))
    result <- run(plain, debug = FALSE, is_plot = FALSE)
    plotted <- run(report, debug = TRUE, is_plot = TRUE)

    expect_equal(plotted, result)
    expect_equal(readLines(paste0(report, "_fCs.txt")), readLines(paste0(plain, "_fCs.txt")))
    expect_length(readLines(paste0(plain, "_Fstat.txt")), 0)
    expect_false(file.exists(paste0(plain, "_F.pdf")))
    expect_gt(file.info(paste0(report, "_F.pdf"))$size, 0)

    log <- readLines(paste0(report, "_Fstat.txt"))
    markers <- if (reference) {
      c("The adj-diff. autocor", "PMT :", "Ignore changepoints", "Common trend")
    } else {
      c("The nominal level", "Common trend", "stepsize=", "PMF :")
    }
    positions <- vapply(markers, function(x) grep(x, log, fixed = TRUE)[1], integer(1))
    expect_false(anyNA(positions))
    expect_equal(positions, sort(positions))

    run(report, debug = FALSE, is_plot = FALSE)
    expect_length(readLines(paste0(report, "_Fstat.txt")), 0)
    expect_equal(readLines(paste0(report, "_fCs.txt")), readLines(paste0(plain, "_fCs.txt")))
  }
})
