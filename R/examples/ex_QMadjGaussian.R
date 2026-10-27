set.seed(1)
ref <- rnorm(120)
base <- ref + rep(c(2, 0), each = 60) + rnorm(120, sd = 0.2)

# %% 无参考序列
r <- QMadjGaussian(base, Ips = c(60, 120), Mq = 3, Iseg.adj = 2, Nadj = 0, Nt = 12)
tapply(r$PA, rep(1:2, each = 60), mean)

# %% 参考序列
r <- QMadjGaussian.wRef(
    base, Dif = base - ref,
    Ips = c(60, 120), Mq = 3, Iseg.adj = 2, Nadj = 0, Nt = 12, Ny4a = 0
)
tapply(r$PA - ref, rep(1:2, each = 60), mean)
