rm_empty <- function(x) {
    if (is.list(x)) {
        x[sapply(x, length) > 0]
    }
    else {
        x[!is.na(x)]
    }
}

range2 <- function(x, y, na.rm = TRUE) {
    c(min(x, y, na.rm = na.rm), max(x, y, na.rm = na.rm))
}

check_dir <- function(path) {
    for (path_i in path){
        if (!dir.exists(path_i)) {
            dir.create(path_i, recursive = TRUE)
        }
    }
    path
}

listk <- function (...) 
{
    cols <- as.list(substitute(list(...)))[-1]
    vars <- names(cols)
    Id_noname <- if (is.null(vars)) 
        seq_along(cols)
    else which(vars == "")
    if (length(Id_noname) > 0) 
        vars[Id_noname] <- sapply(cols[Id_noname], deparse)
    x <- setNames(list(...), vars)
    return(x)
}

#' @export 
num2date <- function(x) {
    if (is(x, "Date")) return(x)
    # else
    as.character(x) %>% 
        gsub("-", "", .) %>% 
        gsub("0000$", "0101", .) %>% 
        gsub(  "00$",   "01", .) %>% as.Date("%Y%m%d")
}

#' @export
date2num <- function(date) {
    as.character(date) %>% gsub("-", "", .) %>%
      as.numeric()
}

is.Date <- function(x) is(x, "Date")
