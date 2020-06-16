test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

sigma <- matrix(c(4,2,2,3), ncol=2)

test_that("rbivnorm works", {
  set.seed(42)
  z1 <- rbivnorm(10000,mean=c(1,2),v1=sigma[1,1],v2=sigma[2,2],v12=sigma[1,2])
  m1 <- colMeans(z1)
  expect_lt(max(abs(m1-c(1,2))),0.1)
  v <- var(z1)
  expect_lt(max(abs(v-sigma)),0.1)
})

test_that("dbivnorm works", {
  set.seed(42)
  z <- rmvnorm(10,mean=c(0,0),sigma=sigma)
  l0 <- dmvnorm(z,mean=c(0,0),sigma=sigma)
  l1 <- dbivnorm(z,mean=c(0,0),v1=sigma[1,1],v2=sigma[2,2],v12=sigma[1,2])
  expect_equal(l0,l1)
})


