library(SimMixDNA)

test_that("Log-Normal model can be constructed", {

  template <- c(1e3,5e2)
  c2 <- 15.

  size_regression <- read_size_regression(system.file("extdata","GlobalFiler_SizeRegression.csv",package = "SimMixDNA"))

  model_settings <- list(locus_names = "locus1",
                         size_regression = size_regression,
                         detection_threshold = c(locus1 = 50))


  model <- log_normal_model(template = template,
                            c2 = c2, model_settings = model_settings)

  expect_identical(model$size_regression, size_regression)
  expect_identical(model$parameters$template, template)
  expect_identical(model$parameters$c2, c2)
})
