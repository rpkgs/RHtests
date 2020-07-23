#' adjust TP according to station meta info
#' @export
adjust_TP <- function(TP, TP_meta, maxgap = 90) {
  TP$date %<>% as.character() %>% gsub("00$|01$", "01", .) %>% as.Date("%Y%m%d")
  # TP_meta = data.table(date = c("1966-11-01", "1976-07-01", "1980-03-01"))
  TP_meta$date %<>% as.Date()

  # 寻找距离最近的日期
  info <- foreach(i = 1:nrow(TP)) %do% {
    diff = difftime(TP$date[i], TP_meta$date, units = "days") %>% as.numeric()
    I = which.min(abs(diff))
    data.table(date_meta = TP_meta$date[I], diff = diff[I])
  } %>% do.call(rbind, .)
  TP %<>% {cbind(.[, 1:9], info)}

  ## 调整Type-1突变点位置
  ## Type-0中仅挑选有metedata支持的TP
  TP0 <- TP[kind == 0 & abs(diff) <= maxgap, ]
  TP0[, date := date_meta]

  TP1 <- TP[kind == 1, ]
  TP1[abs(diff) <= maxgap, `:=`(kind = 0, date = date_meta)]

  TP_final = rbind(TP0, TP1) %>% .[order(date), ]
  TP_final[, date := format(date, "%Y%m00") %>% as.integer()]
  # print(TP_final)
  TP_final
#   if (!all.equal(TP[, 3:9], TP_final[, 3:9])) {
#     StepSize(infile, InCs = "OUTPUT/example01_1Cs.txt", output = "OUTPUT/example01")
#   }
}

#' @param r object returned by [StepSize()]
#' @rdname adjust_TP
#' @export
adjust_step_TP <- function(r) {
  TP2   <- r$turningPoint
  I_del <- TP2[, which.min(abs(stepsize))]
  kind  <- TP2$kind[I_del]
  if (kind == 0) {
    is_keep = TP2[I_del, prob >= probU]
  } else if (kind == 1) {
    is_keep = TP2[I_del, PFx >= PFx95h]
  }
  if (!is_keep) TP2 = TP2[-I_del, ]
  TP2
}
