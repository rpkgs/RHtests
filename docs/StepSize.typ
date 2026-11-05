#set document(title: "StepSize：突变幅度估计与序列订正")
#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "", header: "")

= StepSize：突变幅度估计与序列订正

`StepSize()` 根据给定的候选突变点，估计各点的跳变幅度、复检其显著性，并将序列订正至指定基准段。它不搜索或移动突变点，也不自行删除检验不显著的点；逐步筛选由上层的 `RHtests_stepsize()` 完成。

== 1 输入与分段

输入包括基准序列 `InSeries`、候选突变点表 `InCs` 和输出路径前缀 `output`。`InCs` 至少包含日期 `date` 和类型 `kind`；候选日期映射至该日或之前最近的有效观测，并须按时间排序。突变点属于前一段，下一段从其后一个有效观测开始。

基准段由 `Iadj` 指定：默认 `10000` 表示最后一段，`0` 表示最长段，其他有效正整数表示对应段。已有读取状态时可传入 `.state`，避免重复读取。

== 2 拟合与订正

无参考序列时，处理流程如下：

1. *去季节：*按月日计算平均季节周期并从序列中减去；分别拟合无突变和含给定突变点的趋势。
2. *估计跳变：*含突变点的模型采用共同线性趋势和分段截距；估计趋势时考虑残差的一阶自相关，并迭代更新季节周期。第 $i$ 个突变点的幅度为 $Delta_i = mu_(i+1) - mu_i$，即相邻段截距之差。
3. *均值订正：*对第 $i$ 段加上 $mu_"ref" - mu_i$，使各段对齐基准段，同时保留趋势和季节周期。
4. *分位数订正：*先去除趋势，再调用 `QMadjGaussian()`，最后加回趋势。分位数订正允许不同数值水平采用不同订正量，详见 `QMadj.typ`。
5. *复检：*比较保留和移除各候选点后的拟合残差，计算 PMF 统计量；结合自相关及其上下界对应的概率和临界值，给出显著性判定。

`Mq` 控制分位数组数上限，默认 10；`Mq = 1` 时为单组均值订正，无参考序列的 `QM_adjusted` 列此时不输出该结果。`Ny4a` 控制突变点两侧估计订正量的窗口年数：0 使用整段，正数不超过 5 时按 5 年处理。

=== 2.1 多突变点拟合：`LSmultiple()`

设序列长度为 $n$，`Ips` 为 $(p_1, dots, p_K, p_(K+1))$，其中 $p_(K+1) = n$，$K$ 为突变点数。模型为

$ Y_t = alpha + beta (X_t - overline(X))
  + sum_(j=1)^K gamma_j I(p_j < t <= p_(j+1)) + epsilon_t, quad t = 1, dots, n. $

其中 $Y_t$ 为观测值，$X_t$ 为自变量，$I(A)$ 在条件 $A$ 成立时为 1，否则为 0。各段共用斜率 $beta$；$gamma_j$ 是第 $j+1$ 段相对第一段的截距差，*不是*所有突变点的跳变量。第一处跳变量为 $gamma_1$，此后第 $j$ 处为 $gamma_j - gamma_(j-1)$。$K = 0$ 时求和项为零。

=== 2.2 单突变点拟合：`LSmatrix()`

当 `Ic` 为突变位置 $c$ 时，模型为

$ Y_t = alpha + beta X_t + gamma I(t > c) + epsilon_t, quad t = 1, dots, n. $

其中 $gamma$ 是突变前后的截距差，突变位置 $c$ 属于前一段；`Ic = NA` 时不包含突变项。`LSmatrix()` 直接使用传入的 $X$，不自行中心化；调用方可事先将其减去均值。

== 3 突变点复检

无参考序列使用 PMF 检验。对每个候选点，在其相邻两段的范围内，比较保留与移除该点时的残差平方和。原序列上的 $F$ 统计量乘以位置惩罚因子，得到 `PFx`；对预白化序列分别计算 `prob` 及其上下界 `probL`、`probU`。若预白化计算的 $F$ 为负，代码将 `PFx` 和 `prob` 置为 0。

`prob` 是 $F$ 分布的累积分布函数值，不是常规的上尾 $p$ 值；数值越大对应的统计量越大。`PFx95` 及其上下界来自临界值表，表由指定的置信水平 `plev` 生成，默认 0.95。

#table(
  columns: (auto, 1fr, auto), align: (left, left, left), inset: 5pt,
  [*类型*], [*判定条件，按顺序*], [*Idc*],
  [`kind = 0`], [`PFx >= PFx95h`], [`Yes`],
  [], [`probU < plev`], [`No`],
  [], [`probL < plev`], [`?`],
  [], [其余情况], [`YifD`],
  [`kind = 1`], [`PFx < PFx95l`], [`No`],
  [], [`PFx < PFx95h`], [`?`],
  [], [其余情况], [`Yes`],
)

`YifD` 表示须有可靠元数据支持才能认定显著。`RHtests_stepsize()` 在候选点不少于两个时，每轮仅检查幅度绝对值最小的点：`kind = 0` 时比较 `prob >= probU`，`kind = 1` 时比较 `PFx >= PFx95h`；不满足则删除该点并重新运行，满足则停止。此规则不等同于逐个删除所有 `Idc = No` 的点；只剩一个点时不再筛选。

== 4 有参考序列时

`StepSize.wRef()` 在基准与参考序列的差值上采用 PMT 检验，并将差值分段均值用于相对参考序列的均值订正；同时也根据基准序列自身的分段趋势拟合计算另一套均值订正。

返回的 `stepsize` 来自基准序列自身的分段回归系数，不是差值序列的分段均值之差。分位数订正调用 `QMadjGaussian.wRef()`：按基准值分组，以基准与参考序列的差值估计订正量，并对均值施加约束。

== 5 输出

非 GUI 模式返回 `list(data, TP)`，分别为逐时次序列与突变点表；GUI 模式返回状态码 0。

- `data`：包括原序列 `base`、去季节序列 `base_anomaly`、均值订正结果、趋势与季节拟合结果，以及 `QM_adjusted`。有参考序列时，`base_mean_adj.ref` 按差值估计订正，`base_mean_adj.base` 按基准序列拟合估计订正。
- `TP`：包括突变日期 `date`、类型 `kind`、判定 `Idc`、幅度 `stepsize`、概率和临界值。有参考序列时，`PFx` 等列名沿用无参考版本，但其值对应 PMT 检验。
- `output_fCs.txt`：突变点结果表；`output_Fstat.txt`：拟合和检验摘要，`debug = FALSE` 时为空；`output_F.pdf`：`is_plot = TRUE` 时绘制的诊断图。

#text(size: 9pt)[实现：`R/StepSize.R`、`R/StepSize.wRef.R`；拟合表、检验和日志分别见 `R/stepsize_output.R`、`R/stepsize_test.R` 和 `R/write_log.R`。]
