test_that("multiplication works", {
    infile = system.file("extdata/Example1.dat", package = "RHtests")
    d = Read(infile, "-999.99")
    expect_equal(633, nrow(d))
    expect_true(data.table::is.data.table(d))
})
