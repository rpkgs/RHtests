# %% 
# 两个突变点，三段数据；第 2 段为基准
before <- matrix(c(1, 3, 0), ncol = 1)
after <- matrix(c(2, 5, 0), ncol = 1)
RHtests:::.qm_adjustments(before, after, ref_segment = 2)
# 订正量：第 1 段 +1，第 2 段 0，第 3 段 -2

# %% 
# 两个分位组的订正量分别为 1 和 3
knots <- c(-0.25, 0.25, 0.75, 1.25)
values <- c(1, 1, 3, 3) # 两端延拓
RHtests:::.qm_interpolate(knots, values, prob = c(0.25, 0.5, 0.75))
# 订正量：1、2、3
