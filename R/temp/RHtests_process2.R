RHtests_process2 <- function(data, data.ref = NULL, metadata, prefix = "./OUTPUT/example02",
    maxgap = 90,
    is_plot = TRUE, verbose = TRUE)
{
    check_dir(dirname(prefix))

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
    r$TP %<>% merge_metainfo(metadata)
    r
}
