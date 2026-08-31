test_that("4x fixed-point solver matches closed form", {
  p <- 0.3
  a <- 1/6
  q1 <- get_equilibrium_Q_one_eta(4, p, a)
  q2 <- get_equilibrium_Q_4x_closed(p, a)
  expect_lt(max(abs(q1 - q2)), 1e-10)
})
