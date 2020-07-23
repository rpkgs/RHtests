library(foreach)
library(iterators)
load_all()
Bfile <- system.file("extdata/Example2.dat", package = "RHtests")
Rfile <- system.file("extdata/Example2_Ref.dat", package = "RHtests")

# profvis::profvis({
{
  output  =  "OUTPUT/example02/example02"
  check_dir(output)

  d       <- Read(Bfile)
  TP_meta <- data.table(date = c("1966-11-01", "1976-07-01", "1980-03-01"))
  # profvis::profvis({
    U   <- FindU(NULL, is_plot = FALSE)
    UD  <- FindUD(NULL, InCs = U$turningPoint, output, is_plot = FALSE)

    TP  <- UD$turningPoint
    TP2 <- adjust_TP(TP, TP_meta)
    r   <- StepSize(NULL, TP2, output)
    times <- 1
    while(times < nrow(TP)) {
      TP2 <- adjust_step_TP(r)
      print(TP2)

      if (nrow(TP2) < nrow(r$turningPoint)) {
        times <- times + 1
        # print(times)
        r <- StepSize(infile, TP2, output)
      } else break
    }
}
# })
