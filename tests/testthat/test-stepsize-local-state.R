test_that("direct StepSize calls do not read or overwrite global state", {
  names <- c(
    "Y0", "Nt", "Icy", "Imd", "IY0", "IY0flg", "Ti", "TiB", "IYBflg",
    "itable", "ori.itable", "bdata", "ori.bdata", "ooflg", "owflg", "phi", "PFmax", "PTmax",
    "Nmin", "EB", "ErrorMSG"
  )
  existing_names <- intersect(names, ls(.GlobalEnv, all.names = TRUE))
  previous <- mget(existing_names, envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    rm(list = intersect(names, ls(.GlobalEnv, all.names = TRUE)), envir = .GlobalEnv)
    list2env(previous, envir = .GlobalEnv)
  })
  for (name in names) assign(name, "sentinel", envir = .GlobalEnv)

  points <- data.table::data.table(date = 19740201, kind = 1)
  output <- tempfile()
  on.exit(unlink(Sys.glob(paste0(output, "*"))), add = TRUE)
  base <- system.file("extdata/Example1.dat", package = "RHtests")
  result_noref <- StepSize(base, points, output, is_plot = TRUE, debug = TRUE)
  expect_equal(nrow(result_noref$data), 636)
  expect_equal(nrow(result_noref$TP), 1)
  expect_gt(file.info(paste0(output, "_F.pdf"))$size, 0)

  base <- system.file("extdata/Example2.dat", package = "RHtests")
  reference <- system.file("extdata/Example2_Ref.dat", package = "RHtests")
  result_ref <- StepSize.wRef(base, reference, points, output, is_plot = TRUE, debug = TRUE)
  expect_equal(nrow(result_ref$data), 612)
  expect_equal(nrow(result_ref$TP), 1)
  expect_gt(file.info(paste0(output, "_F.pdf"))$size, 0)

  iterative_noref <- RHtests_stepsize(
    system.file("extdata/Example1.dat", package = "RHtests"), TP2 = points,
    prefix = output, is_plot = FALSE, verbose = FALSE
  )
  iterative_ref <- RHtests_stepsize(
    base, reference, TP2 = points, prefix = output, is_plot = FALSE, verbose = FALSE
  )
  expect_equal(iterative_noref$data, result_noref$data)
  expect_equal(iterative_ref$data, result_ref$data)

  outdir <- tempfile()
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  metadata <- data.table::data.table(date = c("1966-11-01", "1976-07-01", "1980-03-01"))
  full_noref <- RHtests_process(
    system.file("extdata/Example1.dat", package = "RHtests"), metadata = metadata,
    outdir = outdir, is_plot = TRUE, verbose = FALSE
  )
  metadata <- data.table::data.table(date = c("1974-02-01", "1975-11-01"))
  full_ref <- RHtests_process(base, reference, metadata, outdir = outdir, is_plot = TRUE, verbose = FALSE)
  expect_equal(nrow(full_noref$data), 636)
  expect_equal(nrow(full_ref$data), 612)
  daily <- system.file("extdata/RHtests_dlyPrcp_ExampleData.txt", package = "RHtests")
  daily_series <- RHtests:::ReadDLY.g(daily, "-999.99")
  expect_equal(length(daily_series$Y0), length(daily_series$IY0))
  expect_identical(mget(names, envir = .GlobalEnv), setNames(rep(list("sentinel"), length(names)), names))

  # 显式传入状态与直接读取应产生相同的计算结果。
  state <- Read(system.file("extdata/Example1.dat", package = "RHtests"))
  expect_equal(StepSize(InCs = points, output = output, is_plot = FALSE, .state = state), result_noref)
  state <- Read.wRef(base, reference)
  expect_equal(StepSize.wRef(InCs = points, output = output, is_plot = FALSE, .state = state), result_ref)

  empty_points <- data.table::data.table(date = numeric(), kind = numeric())
  result <- StepSize.wRef(base, reference, empty_points, output, is_plot = FALSE, debug = FALSE)
  expect_equal(nrow(result$data), 612)
  expect_length(result$TP, 0)
})

test_that("input errors are reported directly", {
  input <- data.table::fread(system.file("extdata/Example1.dat", package = "RHtests"))
  input[c(1, 2), ] <- input[c(2, 1), ]
  expect_error(Read(input), "not continuous")

  base <- data.table::fread(system.file("extdata/Example2.dat", package = "RHtests"))
  reference <- data.table::fread(system.file("extdata/Example2_Ref.dat", package = "RHtests"))
  reference[3, V3 := 2]
  expect_error(Read.wRef(base, reference), "ref series not continuous")
})
