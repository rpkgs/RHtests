#' Colored terminal output
#' 
#' @param ... Strings to style.
#' @importFrom crayon green red bold underline
#' 
#' @keywords internal
#' @rdname ok
#' @export
ok    <- function(...) cat(green(...))

#' @rdname ok
#' @export
warn  <- function(...) cat(red(...))
# warn  <- function(...) cat(red $ underline (...))

width_str <- function(str, width = NULL){
    if (!is.null(width) && width > 0) {
        pattern <- sprintf("%%%ds", width)
        sprintf(pattern, str)
    } else {
        sprintf("%s", str)
    }
}

num_bad <- function(str, width = NULL, ...){
    str <- width_str(str, width)
    bold $ underline $ red (str)
}

num_good <- function(str, width = NULL, ...){
    str <- width_str(str, width)
    bold $ underline $ green (str)
}
