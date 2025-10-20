context("ridgereg")

test_that("ridgereg returns expected coefficients", {
  skip_if_not_installed("MASS")
  data(mtcars)
  lambda <- 2.5

  fit_custom <- ridgereg(mpg ~ cyl + disp + hp, data = mtcars, lambda = lambda)
  fit_mass <- MASS::lm.ridge(mpg ~ cyl + disp + hp, data = mtcars, lambda = lambda)
  intercept <- fit_mass$ym - sum(fit_mass$xm * fit_mass$coef / fit_mass$scales)
  coef_mass <- c(`(Intercept)` = intercept, fit_mass$coef / fit_mass$scales)

  expect_s3_class(fit_custom, "ridgereg")
  expect_named(coef(fit_custom))
  expect_equal(unname(coef(fit_custom)), unname(coef_mass), tolerance = 1e-2)
})

test_that("predict and fitted values align", {
  data(mtcars)
  model <- ridgereg(mpg ~ cyl + wt + hp, data = mtcars, lambda = 1)

  expect_equal(predict(model)[1:6], fitted(model)[1:6])

  preds_new <- predict(model, newdata = mtcars[1:3, ])
  manual <- as.numeric(model.matrix(~ cyl + wt + hp, mtcars[1:3, ]) %*% coef(model))
  expect_equal(preds_new, manual)
})

test_that("summary returns structured output", {
  data(mtcars)
  model <- ridgereg(mpg ~ cyl + wt + hp, data = mtcars, lambda = 0.5)

  summ <- summary(model)
  expect_s3_class(summ, "summary.ridgereg")
  expect_true(all(c("Estimate", "Std. Error") %in% colnames(summ$coefficients)))
  expect_output(print(summ), "Lambda:")
})

test_that("visualize_airport_delays returns a ggplot", {
  skip_if_not_installed("nycflights13")
  plt <- visualize_airport_delays()
  expect_s3_class(plt, "ggplot")
})
