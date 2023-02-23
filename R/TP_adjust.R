#' @export
merge_metainfo <- function(TP, metadata) {
  if (is.null(TP)) return(NULL)

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
#' @param TP A data.frame, with the columns of `c("kind", "diff", "")`
#' 
#' @return An adjusted date
#' 
#' @export
TP_adjustByMeta <- function(TP, metadata, maxgap = 366) {
  if (is.null(TP)) return(NULL)

  TP %<>% merge_metainfo(metadata)
  
  # Type-0突变点，根据meta信息进行调整
  # (1) 仅挑选有meta支持的TP
  # (2) TP日期根据meta进行微调
  TP0 <- TP[kind == 0 & abs(diff) <= maxgap, ] %>% mutate(date = date_meta)
  
  ## 此举TP1也被干掉，有可能转变为TP0, good option
  TP1 <- TP[kind == 1, ]
  TP1[abs(diff) <= maxgap, `:=`(kind = 0, date = date_meta)]

  if (nrow(TP0) == 0) TP0 = NULL
  if (nrow(TP1) == 0) TP1 = NULL
  TP_final = rbind(TP0, TP1)
  
  ## MERGE DUPLICATED DATE
  TP_final[, .SD[which.min(abs(diff)), ], .(date)] %>% .[order(date), ]
}

#' @rdname TP_adjustByMeta
#' @export
TP_remove_least <- function(TP) {
  # 每次删除掉一个最不显著的突变点
  if (is.null(TP)) return(NULL)
  
  I_del <- TP[, which.min(abs(stepsize))]
  kind  <- TP$kind[I_del]

  if (kind == 0) {
    is_keep = TP[I_del, prob >= probU]
  } else if (kind == 1) {
    is_keep = TP[I_del, PFx >= PFx95h]
  }
  if (!is_keep) TP = TP[-I_del, ]
  TP
}
