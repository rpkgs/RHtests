## 参考站挑选规则

**1. 确定潜在参考站**

Q1: PMF检测无突变点(410 sites in total)

Q2: PMF检测有突变点，但站点位置未发生变迁，采用QM-adjusted结果作为reference

**2. 挑选每个站点的参考站**

1. 350km范围内, corr(diff(Ti), diff(Rj)) > 0.7，时间序列最长的站点

   同时需要满足：|HR -HT| ≤200m if HT<2500 m; |HR -HT| ≤500m if HT≥2500 m

<!-- 2. 如果没有满足的potential reference sites，则采用Q2的站点 -->
2. 如果经过步骤1，仍均无法找到参考站点，则采用无参考站点的矫正方法(PMF test)

**3. PMT实战**

实战

测试reference site出现空值对QM_adj的影响


**城市化站点挑选**

-   突变点显著，RHtests则矫正到Rural (or Urban)，则最终为Rural (or Urban)站点
-   突变点不显著，则采用动态Land Cover划分Urban or Rural
