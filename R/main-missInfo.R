#' merge txt files into a csv
#' 
#' @details 
#' ## WIN errors fixed:
#' [826] 69473, 59034 2429 10839   1851 2019 10 235000002     64      5    106      5 9 9 9 9 9
#'   win_avg: 5000002 -> 20
#' [827] 17129: 53955 3531 11028   4616 2019 11 292500012     33      3     64      3 9 9 9 9 9
#'   win_avg: 500012 -> 12
#' [829]  3792: 51238 4454  8204   5322 2020  1  72500001     21     13     32     13 9 9 9 9 9
#'   win_avg: 2500001 -> 10
#' @param dir_root e.g. `dir_root = "N:/DATA/China/2400climate data" %>% path.mnt()`
#' @export
#' 
#' @examples 
#' \dontrun{
#' dir_root = "N:/DATA/China/2400climate data" %>% path.mnt()
#' merge_mete2000_txts(dir_root, is_save = TRUE)
#' }
merge_mete2000_txts <- function(dir_root, is_save = TRUE) {
    # dir_root = "N:/DATA/China/2400climate data" %>% path.mnt()
    varnames = c("EVP", "GST", "PRE", "PRS", "RHU", "SSD", "TEM", "WIN")
    foreach(varname = varnames[3], i = icount()) %do% {
        outfile = glue("{dir_root}/SURF_CLI_CHN_MUL_DAY_{varname} (195101-202003).csv")
        # if (file.exists(outfile)) return()

        indir = glue("{dir_root}/{varname}")
        files = dir(indir, "*.TXT", full.names = TRUE)

        lst <- foreach(file = files[1:length(files)], i = icount()) %do% {
            runningId(i, 10)
            tryCatch({
                fread(file)
            }, warning = function(e) {
                message(sprintf('[i] %s: %s', i, basename(file), e$message))
            })
        }
        # lst <- llply(files, fread, .progress = "text")
        df = do.call(rbind, lst)
        if (is_save) fwrite(df, outfile)
        invisible()
    }
}
