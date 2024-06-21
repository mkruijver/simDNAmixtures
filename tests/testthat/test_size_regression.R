test_that("Size regression works with AMEL exceptions", {
  filename <- system.file("extdata",
                          "GlobalFiler_SizeRegression.csv",
                          package = "simDNAmixtures")

  regression <- read_size_regression(filename, exceptions = list(
    AMEL = setNames(c(98.5, 104.5), nm = c("X", "Y"))))

  # obtain size for the 12 allele at the vWA locus
  expect_equal(regression("vWA", 12), 112.315164835165 + 12 * 4.0372967032967)

  # verify that AMEL is supported
  expect_equal(regression("AMEL", "X"), 98.5)
  expect_equal(regression("AMEL", "Y"), 104.5)
}
)
