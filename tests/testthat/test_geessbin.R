library(testthat)
library(geessbin)

res <- geessbin(formula = Wheeze ~ City + Time, id = ID,
                repeated = Time, corstr = "ar1", data = wheeze,
                beta.method = "PGEE", SE.method = "MB")
res_sum <- summary(res)

coef1 <- res$coefficients[2]
names(coef1) <- NULL

se1 <- res_sum$coefficients[2, 2]

test_that("Test geessbin not cleared", {
  expect_equal(round(coef1, 3), 0.228)
  expect_equal(round(se1, 4), 0.7498)
})
