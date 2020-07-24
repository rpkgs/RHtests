library(foreach)
library(iterators)
load_all()
Bseries <- system.file("extdata/Example2.dat", package = "RHtests")
Rseries <- system.file("extdata/Example2_Ref.dat", package = "RHtests")

# profvis::profvis({
{
  output  =  "OUTPUT/example02/example02"
  check_dir(output)

  U  <- FindU.wRef(Bseries, Rseries, output)
  UD <- FindUD.wRef(Bseries, Rseries, U$turningPoint, output)

#   d       <- Read(Bfile)
  TP_meta <- data.table(date = c("1974-02-01", "1975-11-01"))
  # profvis::profvis({
    # U   <- FindU(Bfile, Rfile, is_plot = FALSE)
    # UD  <- FindUD(NULL, InCs = U$turningPoint, output, is_plot = FALSE)
    TP  <- UD$turningPoint
    TP2 <- adjust_TP(TP, TP_meta, maxgap = 90)
    r   <- StepSize.wRef(Bseries, Rseries, TP2, output)
    times <- 1
    while(times < nrow(TP)) {
      TP2 <- adjust_step_TP(r)

      if (nrow(TP2) < nrow(r$turningPoint)) {
        print(TP2)
        times <- times + 1
        # print(times)
        r <- StepSize.wRef(Bseries, Rseries, TP2, output)
      } else break
    }
}
# })
