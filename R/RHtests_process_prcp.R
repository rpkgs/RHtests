# RHtests_stepsize_prcp(data, TP2, prefix)
#' @export
RHtests_stepsize_prcp <- function(data = NULL, TP2,
    prefix = "./OUTPUT/example03", is_plot = FALSE, verbose = TRUE)
{
  # RHtests_read(data, data.ref)
  # has_ref = !is.null(data.ref)
  # FUN_step <- if (has_ref) StepSize.wRef else StepSize

  r = StepSize.dlyPrcp(data, TP2, prefix)
  # r   <- FUN_step(InCs = TP2, output = prefix, is_plot = is_plot)

  times <- 1
  while (times < nrow(TP2)) {
    if (length(r$TP) == 0) return(NULL)
    TP2 <- adjust_step_TP(r)

    if (nrow(TP2) < nrow(r$TP)) {
      if (verbose) print(TP2)
      times <- times + 1
      r = StepSize.dlyPrcp(data, TP2, prefix)
      # r <- FUN_step(InCs = TP2, output = prefix, is_plot = is_plot)
    } else {
      break
    }
  }
  r
}

#' @export
RHtests_process_prcp <- function(data, metadata, prefix = "./OUTPUT/example02",
    maxgap = 366,
    is_plot = TRUE, verbose = TRUE)
{
  check_dir(dirname(prefix))

  # has_ref = !is.null(data.ref)
  # FUN_FindU  <- if (has_ref) FindU.wRef else FindU
  # FUN_FindUD <- if (has_ref) FindUD.wRef else FindUD

  # RHtests_read(data, data.ref)
  U <- FindU.dlyPrcp(data, prefix)
  if (is_empty(U$TP)) return(NULL)

  UD <- FindUD.dlyPrcp(data, U$TP, prefix)
  if (is_empty(UD$TP)) return(NULL)

  TP  <- UD$TP
  TP2 <- adjust_TP(TP, metadata, maxgap = maxgap)

  # StepSize.dlyPrcp(data, TP2, prefix)
  # r <- RHtests_stepsize(data = NULL, data.ref = NULL, TP2, has_ref,
  #     prefix, is_plot, verbose)
  r <- RHtests_stepsize_prcp(data, TP2, prefix)
  r$TP %<>% merge_metainfo(metadata)
  r
}
