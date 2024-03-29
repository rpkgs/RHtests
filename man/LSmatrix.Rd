% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LSmultiRed.R
\name{LSmatrix}
\alias{LSmatrix}
\title{single-breakpoint piecewise regression (same slope, different interception)}
\usage{
LSmatrix(Y, T_anorm, Ic, ...)
}
\arguments{
\item{Y}{The response vector}

\item{Ic}{the position of change point, \verb{[1, Ic]} and \verb{[Ic + 1, end]} are
corresponding to the perid

#' @return
\itemize{
\item \code{sig}   : coefficients
\item \code{fitted}: fitted value
\item \code{resi}  : residual of fitted value
\item \code{SSE}   : sum of square error
}}

\item{T}{The predictor vector}
}
\description{
single-breakpoint piecewise regression (same slope, different interception)
}
