#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

#**** Copyright agreement ****
#       All users of any code in this package agree to the terms and conditions described in the 
#file Copyright_RHtests__RClimDex_SoftwarePackages.pdf, which is also included in this package.
#*****************************
# last updated at 2019-03-01
# in QMadjGaussian.wRef(), change osmean estimation from the way QMadjGaussian to just take the
# mean value of adjusted segments, since there has extra adjustments for the difference b/w
# QMadj and mean-adj

# last updated at 2018-03-01
# in QMadjGaussian.wRef(), added check sample size during EBa calculation, and reduce Mq once
#   sample size too small detected (the same consequence as EBb estimation)

# last updated at 2016-12-28
# set all QMadj output as missing once Mq==1

# last updated at 2016-10-05
# changed FindU.wRef(), FindUD.wRef(), StepSize.wRef(), switched LSmultipleRed to LSmultiRedCycle
# for trend estimation on original and mean-adjusted series

# last updated at 2014-02-03
# 1) Fixed Read.wRef(), changed merge() parameter sort=F to sort=T for itable.nmb, otherwise, 
# non-match date for Bseries will be at the end of itable.nmb
# 2) Changed all pdf(), added paper='letter'
# 3) Changed OnQMadjGaussian.wRef(), make output filename correctly involved reference file name.

# last updated on 2013-06-28
# Changed all as.real() to as.numeric(), since R 2.15 no longer support as.real()

# last updated on 2013-04-03
# Changed Mq upper limit in StartGUI(), QMadjGaussian(), QMadjGaussian.wRef(), from 20 to 100,
# changed button locations on main manu.

# last updated on 2013-04-03
# Changed plots for QMadjust distribution

# last updated on 2013-03-19
# Added function QMadjGaussian.wRef(), changed FindU.wRef(), FindUD.wRef(), 
#   StepSize.wRef(), replaced QMadjGaussian() with QMadjGaussian.wRef().
#   Modified Read.wRef(), append ref series to bdata matrix.
# Added a button 'QMadj.wRef' at StartGUI(), corresponding subroutine: QMadj.GaussianDLY.wRef

# last updated on 2012-09-25
# fixed a bug at StepSize() and StepSize.wRef() regarding to fitted value for Feb. 29th 

# last updated on 2012-04-26
# at FindU() FindUD() StepSize() FindU.wRef() FindUD.wRef() and 
# StepSize.wRef(), add fitted value for Feb. 29th to fitted 
# value data files.

# last updated on 2011-12-11
# Change output format for *Ustat.txt in *.wRef(), from "PMF" -> "PMT"
# Change QMadj.GuassianDLY() and ReadDLY.g(), add vector olflg which indicate
#   non-missing data only, but contains Feb 29th data.

# last updated on 2010-06-08
#   In FindU(), FindUD(), StepSize() and 3 corresponding .wRef() functions,
#   added a text indication "Yes ", "YifD", "No  " or "?   " in changepoint list
#   output file, also changed corresponding read-in part.

# last updated on 2010-10-29
# add a sample format for 1Cs when there has no changepoint
# last updated on 2010-05-06
# in QMadjGaussian(), add Nadj option, set empirical prob. segment as length
# Nadj instead of whole segment.
#
# Bug fixed at 2010-02-08
# In FindU.wRef(), FindUD.wRef(), if 0 changepoint found, set meanhatD with
#   fitted value, but not 0. Output *_U.dat or *_UD.dat affected.

# NT=1, 12, and 365 for annual, monthly, and daily series, respectively
# p.lev is the nominal level of confidence; (1 - p.lev) is the nominal level
#     of significance
# choose one of the following 6 p.lev values: 0.75, 0.80, 0.90, 0.95, 0.99,
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
