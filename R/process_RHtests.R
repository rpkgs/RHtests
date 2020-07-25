#' process RHtests without reference series
#'
#' @param data A data.frame or data.table, with the columns of c('year', 'month', 'day', 'data')
#' @param data.ref A data.frame or data.table, with the columns of c('year', 'month', 'day', 'data').
#' The reference series of `data`, should have the same length as `data`.
#' @param metadata A data.frame or data.table, with the columns of `TurningPoint` date.
#'
#' @example man/examples/run_ex01.R
#' @example man/examples/run_ex02.R
#'
#' @export
process_RHtests <- function(data, data.ref = NULL, metadata, prefix = "./OUTPUT/example02",
    maxgap = 90,
    is_plot = TRUE, verbose = TRUE)
{
    check_dir(dirname(prefix))

    has_ref = !is.null(data.ref)
    FUN_FindUD <- if (has_ref) FindUD.wRef else FindUD
    FUN_step <- if (has_ref) StepSize.wRef else StepSize

    if (has_ref) {
        Read.wRef(data, data.ref)
        U  <- FindU.wRef(NULL, NULL, prefix, is_plot = is_plot)        
    } else {
        d <- Read(data)
        U <- FindU(output = prefix, is_plot = is_plot)
    }
    if (is_empty(U$turningPoint)) return(NULL)

    UD <- FUN_FindUD(InCs = U$turningPoint, output = prefix, is_plot = is_plot)
    if (is_empty(UD$turningPoint)) return(NULL)
    
    TP  <- UD$turningPoint
    TP2 <- adjust_TP(TP, metadata, maxgap = maxgap)
    r   <- FUN_step(InCs = TP2, output = prefix, is_plot = is_plot)

    times <- 1
    while (times < nrow(TP)) {
        if (length(r$turningPoint) == 0) return(NULL)
        TP2 <- adjust_step_TP(r)

        if (nrow(TP2) < nrow(r$turningPoint)) {
            if (verbose) print(TP2)
            times <- times + 1
            r <- FUN_step(InCs = TP2, output = prefix, is_plot = is_plot)
        } else {
            break
        }
    }
    r
}

is_empty <- function(x) length(x) == 0

#' @import ggplot2
#' @export
plot_RHtests <- function(r, outfile = "RHtests-noref.pdf", is_plotly = FALSE) {
    data = r$data
    d2 <- melt(data, c("id", "date", "base"))
    # d2[value < -100, value := NA]
    p <- ggplot(d2, aes(id, value, color = variable)) +
        geom_line(data = data, aes(id, base, color = NULL), color = alpha("black", 0.6)) +
        geom_line()

    latticeGrob::write_fig(p, outfile, 10, 6)
    if (is_plotly) plotly::ggplotly(p)
}
