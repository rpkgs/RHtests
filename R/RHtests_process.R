RHtests_read <- function(data, data.ref = NULL, plev = 0.95) {
  if (is.null(data)) return()
  if (is.null(data.ref)) {
    Read(data, plev = plev)
  } else {
    Read.wRef(data, data.ref, plev = plev)
  }
}


#' RHtests_stepsize
#' stepsize will not alter TP2, only delte non-significant record.
#'
#' @export
RHtests_stepsize <- function(data = NULL, data.ref = NULL, TP2,
    has_ref = !is.null(data.ref),
    prefix = "./OUTPUT/example03", is_plot = FALSE, verbose = TRUE)
{
    mkdir(prefix)
    
    RHtests_read(data, data.ref)
    # has_ref = !is.null(data.ref)
    FUN_step <- if (has_ref) StepSize.wRef else StepSize

    if (nrow(TP2) == 0) return(NULL)
    r   <- FUN_step(InCs = TP2, output = prefix, is_plot = is_plot)

    times <- 1
    while (times < nrow(TP2)) {
        if (!is.list(r) || length(r$TP) == 0) return(NULL)
        TP2 <- adjust_step_TP(r)

        if (nrow(TP2) < nrow(r$TP)) {
            if (verbose) print(TP2)
            times <- times + 1
            r <- FUN_step(InCs = TP2, output = prefix, is_plot = is_plot)
        } else {
            break
        }
    }
    r
}

#' process RHtests without reference series
#'
#' @param data A data.frame or data.table, with the columns of c('year', 'month', 'day', 'data')
#' @param data.ref A data.frame or data.table, with the columns of c('year', 'month', 'day', 'data').
#' The reference series of `data`, should have the same length as `data`.
#' @param metadata A data.frame or data.table, with the columns of `TurningPoint` date.
#'
#' @example R/examples/run_ex01.R
#' @example R/examples/run_ex02.R
#'
#' @export
RHtests_process <- function(data, data.ref = NULL, metadata, prefix = "./OUTPUT/example02",
    maxgap = 366,
    is_plot = TRUE, verbose = TRUE)
{
    mkdir(dirname(prefix))

    has_ref = !is.null(data.ref)
    FUN_FindU  <- if (has_ref) FindU.wRef else FindU
    FUN_FindUD <- if (has_ref) FindUD.wRef else FindUD

    RHtests_read(data, data.ref)
    U  <- FUN_FindU(output = prefix, is_plot = is_plot)
    if (is_empty(U$TP)) return(NULL)

    UD <- FUN_FindUD(InCs = U$TP, output = prefix, is_plot = is_plot)
    if (is_empty(UD$TP)) return(NULL)

    TP  <- UD$TP
    TP2 <- adjust_TP(TP, metadata, maxgap = maxgap)

    r <- RHtests_stepsize(data = NULL, data.ref = NULL, TP2, has_ref,
        prefix, is_plot, verbose)
    # if null returned,
    if (is.null(r$TP) || nrow(r$TP) == 0) return(NULL)
    r$TP %<>% merge_metainfo(metadata)
    r
}

RHtests_process2 <- function(data, data.ref = NULL, metadata, prefix = "./OUTPUT/example02",
    maxgap = 90,
    is_plot = TRUE, verbose = FALSE)
{
    mkdir(dirname(prefix))

    has_ref = !is.null(data.ref)
    FUN_FindU  <- ifelse(has_ref, FindU.wRef, FindU)
    FUN_FindUD <- ifelse(has_ref, FindUD.wRef, FindUD)

    RHtests_read(data, data.ref)
    U  <- FUN_FindU(output = prefix, is_plot = is_plot)
    if (is_empty(U$TP)) return(NULL)

    UD <- FUN_FindUD(InCs = U$TP, output = prefix, is_plot = is_plot)
    if (is_empty(UD$TP)) return(NULL)

    TP  <- UD$TP
    TP2 <- adjust_TP(TP, metadata, maxgap = maxgap)

    r <- RHtests_stepsize(data = NULL, data.ref = NULL, TP2, has_ref,
        prefix, is_plot, verbose)
    r$TP %<>% merge_metainfo(metadata)
    r
}
