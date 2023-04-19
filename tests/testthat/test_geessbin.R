library(testthat)
library(geessbin)

res <- geessbin(formula = Wheeze ~ City + factor(Age), data = wheeze, id = ID,
                corstr = "ar1", repeated = Age, beta.method = "PGEE",
                SE.method = "MB")

res_sum <- summary(res)

coef1 <- res$coefficients[2]
names(coef1) <- NULL

se1 <- res_sum$coefficients[2, 2]

test_that("Test geessbin not cleared", {
  expect_equal(round(coef1, 3), 0.226)
  expect_equal(round(se1, 4), 0.8274)
})
