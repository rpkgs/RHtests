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

