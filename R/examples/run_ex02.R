# %% example 02
read_data <- \(f) {
    data <- data.table::fread(f, na.strings = "-999.99")
    setNames(data, c("year", "month", "day", "data"))
}

base <- system.file("extdata/Example2.dat", package = "RHtests") |> read_data()
ref <- system.file("extdata/Example2_Ref.dat", package = "RHtests") |> read_data()

metadata <- data.table::data.table(date = c("1974-02-01", "1975-11-01"))

r <- RHtests_process(base, ref, metadata, outdir = "OUTPUT/example02")
plot_RHtests(r, "RHtests-ref.pdf")

# expect_equal(nrow(r$turningPoint), 11)
# expect_equal(as.character(r$turningPoint[kind == 0]$date), metadata$date)
# expect_equal(nrow(r$data), 612)
