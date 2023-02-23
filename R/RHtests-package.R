#' @keywords internal
#' @import data.table magrittr
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

# **** Copyright agreement ****
# All users of any code in this package agree to the terms and conditions
# described in the file Copyright_RHtests__RClimDex_SoftwarePackages.pdf, which
# is also included in this package.
# *****************************

# NT=1, 12, and 365 for annual, monthly, and daily series, respectively
# plev is the nominal level of confidence; (1 - plev) is the nominal level
#     of significance
# choose one of the following 6 plev values: 0.75, 0.80, 0.90, 0.95, 0.99,
#     0.9999
# Mq (=0, 1, 2, ..., 20) is the number of points (categories) on which
#     the empirical cumulative distribution function are estimated.
#     If Mq=0, the actual Mq is determined by the length of the shortest segment
#     Mq=1 corresponds to mean adjustments (one adjustment for all data in 
#     the same segment)
# In the output: Ns >= 0 is the number of changepoints identified;
#     changepoint positions are Ip(1), Ip(2),...,Ip(Ns)
# If Iadj = 0, the data series is adjusted to the longest segment;
#     otherwise the data series is adjusted to the chosen segment Iadj (if the 
#     given integer Iadj > Ns+1, we re-set Iadj=Ns+1, which corresponds to
#     adjusting the series to the last segment). Set Iadj = 10000 if you want 
#     to ensure that the series is adjusted to the last segment

#' @importFrom utils write.table read.table
#' @importFrom graphics axis lines par text title
#' @importFrom grDevices dev.off pdf
#' @importFrom stats .lm.fit lm pf pt qt var
.onLoad <- function(libname, pkgname) {
  if (getRversion() >= "2.15.1") {
    utils::globalVariables(
      c(
        ".", ".SD", ".N",
        "ooflg", "owflg", 
        "bdata", "ori.bdata", 
        "itable", "ori.itable", 
        "ErrorMSG", 
        "olflg", "Nmin", 
        "Ti", "Y0", "Nt", "IY0", "IY0flg", "Icy", "Imd", "phi", "PFmax", "PTmax"
      )
    )
  }
}
