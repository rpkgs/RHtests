# 拟合基准序列，交替更新分段均值和季节周期。
fit_reference_cycle <- function(base, seasonal, previous_cycle, cycle_index, cycle_values,
                                time, points, continuous, reference_segment, n_cycles) {
  segments <- segment_bounds(points)
  n_breaks <- length(points) - 1
  iteration <- 0
  repeat {
    iteration <- iteration + 1
    fit <- LSmultipleRed(base, time, points, continuous, n_cycles)
    if (n_breaks == 0) {
      cycle <- seasonal
      break
    }

    # 暂时将各段调整至基准段，以估计共同的季节周期。
    delta <- fit$sig[-1] - fit$sig[reference_segment + 1]
    base <- add_segment_adjustment(base + seasonal, segments, -delta)
    cycle <- map_dbl(cycle_values, \(value) mean(base[cycle_index == value], na.rm = TRUE))
    seasonal <- cycle[match(cycle_index, cycle_values)]
    base <- add_segment_adjustment(base - seasonal, segments, delta)

    variability <- sqrt(var(cycle))
    converged <- is.na(variability) ||
      max(abs(previous_cycle - cycle)) < variability / 1000 || iteration > 20
    previous_cycle <- cycle
    if (converged) break
  }
  list(fit = fit, base = base, seasonal = seasonal, cycle = cycle)
}
