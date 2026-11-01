test_that("multiplication works", {
    infile = system.file("extdata/Example1.dat", package = "RHtests")
    state <- Read(infile, "-999.99")
    expect_true(S7::S7_inherits(state, RHtests:::RHState))
    expect_equal(633, length(state@series$Y0))
    expect_true(S7::S7_inherits(Read(state@series$ori.itable), RHtests:::RHState))

    # dat_example01 <- ori.itable %>% as.data.table()
    # use_data(dat_example01, overwrite = TRUE)
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
