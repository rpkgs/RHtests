test_that("reference candidate search preserves the expected changepoints", {
  base <- system.file("extdata/Example2.dat", package = "RHtests")
  ref <- system.file("extdata/Example2_Ref.dat", package = "RHtests")
  state <- Read.wRef(base, ref)
  Y0 <- state@series$Y0
  cv <- state@critical_values

  result <- search_pmt_candidates(Y0, length(Y0), 0, 5, cv$phi, cv$PTmax, 0.95)
  expect_equal(result$Ips, c(
    74, 139, 160, 166, 177, 230, 242, 265, 271, 276, 285, 291, 301, 314,
    331, 343, 373, 382, 389, 447, 457, 469, 485, 500, 507, 513, 519, 605, 610
  ))
  expect_equal(result$Ns, 28)
})
