# 写出有参考序列时的输入信息及差值序列自相关。
log_ref_inputs <- function(file, base, reference, fit, probability) {
  cat(paste("Input Base Series:", base, "\n"), file = file)
  cat(paste("Input Ref Series:", reference, "\n"), file = file, append = TRUE)
  fprintf("The adj-diff. autocor is: %.4f ( %.4f , %.4f p= %.4f )\n",
    fit$cor, fit$corl, fit$corh, probability, file = file, append = TRUE
  )
}

# 写出有参考序列时订正前后的趋势拟合结果。
log_ref_fits <- function(file, fit_base, fit_adjusted, mean_difference) {
  fprintf(
    "Ignore changepoints -> trend0 = %.6f ( %.6f , %.6f ) (p= %.4f ); cor= %.4f ( %.4f , %.4f )\n\n",
    fit_base$trend, fit_base$betaL, fit_base$betaU, fit_base$p.tr,
    fit_base$cor, fit_base$corl, fit_base$corh, file = file, append = TRUE
  )
  cat("Step-sizes estimated from difference series:\n", file = file, append = TRUE)
  if (length(mean_difference) > 1) {
    cat(round(diff(mean_difference), 4), file = file, append = TRUE, fill = 80)
  }
  fprintf(
    "\n after such adjustments, the base series trend= %.6f ( %.6f , %.6f ) (p= %.4f ); cor= %.3f ( %.3f , %.3f )\n\n",
    fit_adjusted$trend, fit_adjusted$betaL, fit_adjusted$betaU, fit_adjusted$p.tr,
    fit_adjusted$cor, fit_adjusted$corl, fit_adjusted$corh, file = file, append = TRUE
  )
}

# 写出有参考序列时基准序列的趋势和突变幅度。
log_ref_result <- function(file, fit, probability, n_breaks, stepsizes) {
  cat("Common trend TPR fit to the de-seasonalized Base series:\n", file = file, append = TRUE)
  fprintf("#steps=  %d ; trend= %.6f (p= %.4f ); cor= %.4f ( %.4f , %.4f )  p= %.4f \n",
    n_breaks, fit$trend, fit$p.tr, fit$cor, fit$corl, fit$corh, probability,
    file = file, append = TRUE
  )
  if (n_breaks > 0) cat(round(stepsizes, 4), file = file, append = TRUE, fill = 80)
}

# 整理并写出突变点表，保留两种流程原有的文件格式。
write_points_table <- function(points, file, reference = FALSE, stepsizes = NULL) {
  if (!length(points)) return(points)
  if (!reference) fprintf("#  %d changepoints in Series \n", length(points), file = file)
  points <- do.call(rbind, points)
  if (reference) {
    points$stepsize <- stepsizes
    columns <- 4:ncol(points)
  } else {
    columns <- 6:13
  }
  for (column in columns) points[[column]] <- round(points[[column]], 4)
  fwrite(points, file, sep = if (reference) "," else "\t", append = TRUE, col.names = TRUE)
  points
}
# 写出突变点检验结果；PMF 和 PMT 使用各自的统计量名称。
write_test_log <- function(results, method, file) {
  statistic_name <- if (method == "PMF") "PFmax" else "PTmax"
  format <- paste0(
    method, " : c=%4d; (Time %10d); Type= %4d; p=%10.4f",
    "(%10.4f-%10.4f); ", statistic_name,
    "=%10.4f; CV95=%10.4f(%10.4f-%10.4f); Nseg=%4d\n"
  )
  for (result in results) {
    fprintf(format, result$Ic, result$date, result$kind, result$prob,
      result$probL, result$probU, result$PFx, result$PFx95,
      result$PFx95l, result$PFx95h, result$Nseg, file = file, append = TRUE
    )
  }
}

# 写出相对于基准段的分位数匹配诊断结果。
write_qm_log <- function(QMout, segments, Iseg.adj, Ny4a, adjustment, file) {
  fprintf("Nseg_shortest = %d ; Mq =  %d ; Ny4a =  %g \n",
    QMout$Nseg.mn, QMout$Mq, Ny4a, file = file, append = TRUE
  )
  fprintf("\n Adjust to segment %d : from %d to %d \n",
    Iseg.adj, segments$start[Iseg.adj], segments$end[Iseg.adj], file = file, append = TRUE
  )
  if (QMout$Mq <= 1) return(invisible(NULL))

  Ns <- length(segments$end) - 1
  oline <- paste("#Fcat: frequency category boundaries\n",
    "#DP: Difference in the category means\n#", sep = ""
  )
  for (i in seq_len(Ns)) oline <- paste(oline, paste0("Fcat.Seg", i), paste0("DP.Seg", i))
  cat(paste(oline, "\n"), file = file, append = TRUE)
  write.table(round(QMout$osmean, 4), file = file, append = TRUE, row.names = FALSE, col.names = FALSE)
  for (i in seq_along(segments$end)) {
    if (i != Iseg.adj) {
      fprintf("Seg. %d: mean of QM-adjustments =%.4f\n", i, adjustment(i), file = file, append = TRUE)
    }
  }
}
