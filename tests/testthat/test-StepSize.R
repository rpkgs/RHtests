test_that("StepSize with noref works", {
  infile <- system.file("extdata/Example1.dat", package = "RHtests")
  output = "../../OUTPUT/example01/example01"
  check_dir(dirname(output))

  d <- Read(infile)
  TP_meta <- data.table(date = c("1966-11-01", "1976-07-01", "1980-03-01"))
  # profvis::profvis({
  U <- FindU(NULL, output, is_plot = FALSE)
  UD <- FindUD(NULL, InCs = U$turningPoint, output, is_plot = FALSE)

  TP <- UD$turningPoint
  TP2 <- adjust_TP(TP, TP_meta)
  r <- StepSize(NULL, TP2, output)
  times <- 1
  while (times < nrow(TP)) {
    TP2 <- adjust_step_TP(r)
    
    if (nrow(TP2) < nrow(r$turningPoint)) {
      print(TP2)
      times <- times + 1
      # print(times)
      r <- StepSize(infile, TP2, output)
    } else {
      break
    }
  }
  expect_equal(nrow(r$turningPoint), 3)
})
