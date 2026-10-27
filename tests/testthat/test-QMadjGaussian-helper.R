test_that("cumulative adjustments use the chosen reference segment", {
  before <- rbind(c(1, 2), c(3, 5), c(0, 0))
  after <- rbind(c(2, 4), c(5, 8), c(0, 0))
  expected <- rbind(c(1, 2), c(0, 0), c(-2, -3))
  expect_equal(RHtests:::.qm_adjustments(before, after, 2), expected)
  expect_equal(RHtests:::.qm_adjustments(before[, 1, drop = FALSE],
    after[, 1, drop = FALSE], 2), expected[, 1, drop = FALSE])
})

test_that("natural interpolation preserves legacy boundary and extrapolation", {
  for (Mq in c(1, 2, 10)) {
    knots <- (0:(Mq + 1) - .5) / Mq
    values <- sin(seq_along(knots))
    values[c(1, Mq + 2)] <- values[c(2, Mq + 1)]
    prob <- c(knots[1] - 1, knots, 0, .5, 1, knots[Mq + 2] + 1)
    derivatives <- RHtests:::splineN(knots, values, 2E30, 2E30)
    expected <- vapply(prob, function(p) RHtests:::splintN(knots, values, derivatives, p), numeric(1))
    expect_equal(RHtests:::.qm_interpolate(knots, values, prob), expected, tolerance = 1e-12)
  }
})
