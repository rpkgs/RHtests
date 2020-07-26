#' @export
merge_metainfo <- function(TP, metadata) {
  TP$date %<>% num2date()
  if (nrow(metadata) == 0 || length(metadata) == 0) {
    return(cbind(TP[, 1:9], date_meta = NA, diff = 9999))
  }
  metadata$date %<>% num2date()

  # 寻找距离最近的日期
  info <- foreach(i = 1:nrow(TP)) %do% {
    diff = difftime(TP$date[i], metadata$date, units = "days") %>% as.numeric()
    I = which.min(abs(diff))
    data.table(date_meta = metadata$date[I], diff = diff[I])
  } %>% do.call(rbind, .)
  TP %<>% {cbind(.[, 1:9], info)}
}

#' adjust TP according to station meta info
#'
#' @return An adjusted date
#'
#' @export
adjust_TP <- function(TP, metadata, maxgap = 366) {
  TP %<>% merge_metainfo(metadata)
  ## 调整Type-1突变点位置
  ## Type-0中仅挑选有metedata支持的TP
  TP0 <- TP[kind == 0 & abs(diff) <= maxgap, ]
  TP0[, date := date_meta]
  ## 补充metadata 的TP，stepsize为进一步剔除不显著部分

  ## 此举TP1也被干掉，有可能转变为TP0, good option
  TP1 <- TP[kind == 1, ]
  TP1[abs(diff) <= maxgap, `:=`(kind = 0, date = date_meta)]

  if (nrow(TP0) == 0) TP0 = NULL
  if (nrow(TP1) == 0) TP1 = NULL
  TP_final = rbind(TP0, TP1)
  ## MERGE DUPLICATED DATE
  TP_final %<>% plyr::ddply(.(date), function(d) {
    as.data.table(d)[which.min(abs(diff)), ]
  }) %>% data.table()
  TP_final[order(date), ]
}

#' @param r object returned by [StepSize()]
#' @rdname adjust_TP
#' @export
adjust_step_TP <- function(r) {
  TP2   <- r$TP
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
