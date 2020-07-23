library(foreach)
library(iterators)
load_all()
infile <- system.file("extdata/Example1.dat", package = "RHtests")

profvis::profvis({
{
  output  =  "OUTPUT/example01/example01"
  check_dir(output)

  d       <- Read(infile)
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
})

# TP1 %<>% cbind(info)
library(ggplot2)
library(latticeGrob)
d = fread("OUTPUT/example01_F.dat")
d[base < -100, base := NA]
d2 <- d %>% melt(c("id", "date", "base"))
d2[value < -100, value := NA]
p <- ggplot(d2, aes(id, value, color = variable)) +
  geom_line(data = d, aes(id, base, color = NULL), color = alpha("black", 0.6)) +
  geom_line()
write_fig(p, "a.pdf", 10, 6)
plotly::ggplotly(p)
