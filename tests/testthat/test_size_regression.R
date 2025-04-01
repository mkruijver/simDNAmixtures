test_that("Size regression works with AMEL exceptions", {
  filename <- system.file("extdata",
                          "GlobalFiler_SizeRegression.csv",
                          package = "simDNAmixtures")

  regression <- read_size_regression(filename, exceptions = list(
    AMEL = stats::setNames(c(98.5, 104.5), nm = c("X", "Y"))))

  # obtain size for the 12 allele at the vWA locus
  expect_equal(regression("vWA", 12), 112.315164835165 + 12 * 4.0372967032967)

  # verify that AMEL is supported
  expect_equal(regression("AMEL", "X"), 98.5)
  expect_equal(regression("AMEL", "Y"), 104.5)
}
)

test_that("Size regression works with repeat length by marker", {
  filename <- system.file("extdata",
                          "GlobalFiler_SizeRegression.csv",
                          package = "simDNAmixtures")

  gf <- gf_configuration()

  regression <- read_size_regression(filename, repeat_length_by_marker =
                                       gf$repeat_length_by_marker)

  # obtain size for the 11.3 allele at the D2S441 locus
  expect_equal(regression("D2S441", 11.3),
               44.7888995552663 + 11.75 * 4.05775572187873)

  # obtain size for the 12 allele at the D2S441 locus
  expect_equal(regression("D2S441", 12),
               44.7888995552663 + 12 * 4.05775572187873)

}
)
