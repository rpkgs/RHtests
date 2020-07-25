test_that("multiplication works", {
    infile = system.file("extdata/Example1.dat", package = "RHtests")
    d = Read(infile, "-999.99")
    expect_equal(633, nrow(d))
    expect_true(data.table::is.data.table(d))

    expect_true({
      Read(ori.itable)
      TRUE
    })

    # dat_example01 <- ori.itable %>% as.data.table()
    # use_data(dat_example01, overwrite = TRUE)

    ## 2.
    Bseries <- system.file("extdata/Example2.dat", package = "RHtests")
    Rseries <- system.file("extdata/Example2_Ref.dat", package = "RHtests")

    Read.wRef(Bseries, Rseries)

    base <- fread(Bseries, na.strings ="-999.99") %>%
      set_colnames(c("year", "month", "day", "data"))
    ref <- fread(Rseries, na.strings ="-999.99") %>%
      set_colnames(c("year", "month", "day", "data.ref"))
    # dat_example02 <- merge(base, ref)
    # use_data(dat_example02, overwrite = TRUE)
    Read.wRef(base, ref)
})
