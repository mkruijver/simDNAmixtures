
test_that("Wide references (matrix) convert to allele tables", {

  x <- structure(c("18", "17", "17", "17", "15", "17", "18", "15", "16",
                   "18", "18", "18", "16", "18", "16", "16", "10", "9", "9", "11",
                   "13", "11", "13", "10"),
                 dim = c(4L, 6L), dimnames = list(c("1","2", "S1", "S2"),
                                                  c("D3S1358", "D3S1358(2)", "vWA", "vWA(2)",
                                                    "D16S539", "D16S539(2)")))

  # convert
  tabs <- .wide_references_to_allele_tables(x)

  expect_true(is.character(tabs[[1]]$Allele1))

  # verify dimensions
  number_of_samples_in <- nrow(x)
  number_of_samples_out <- length(tabs)
  expect_identical(number_of_samples_in, number_of_samples_out)

  number_of_loci_in <- ncol(x) / 2
  number_of_loci_in_per_sample <- rep(number_of_loci_in, number_of_samples_in)
  number_of_loci_out_per_sample <- unname(sapply(tabs, nrow))
  expect_equal(number_of_loci_in_per_sample, number_of_loci_out_per_sample)

  # verify a couple of datapoints
  expect_identical(x["S1", "vWA"], tabs$S1$Allele1[2])
  expect_identical(x["S2", "vWA(2)"], tabs$S2$Allele2[2])
})


test_that("Wide references (df) convert to allele tables", {

  wide_in <- structure(list(D3S1358 = c("18", "17", "17", "17"),
                            `D3S1358(2)` = c("15","17", "18", "15"),
                            vWA = c("16", "18", "18", "18"),
                            `vWA(2)` = c("16","18", "16", "16"),
                            D16S539 = c("10", "9", "9", "11"),
                            `D16S539(2)` = c("13","11", "13", "10")),
                       class = "data.frame",
                       row.names = c("1", "2", "S1", "S2"))

  # convert
  tabs <- .wide_references_to_allele_tables(wide_in)

  expect_true(is.character(tabs[[1]]$Allele1))

  # verify dimensions
  number_of_samples_in <- nrow(wide_in)
  number_of_samples_out <- length(tabs)
  expect_identical(number_of_samples_in, number_of_samples_out)

  number_of_loci_in <- ncol(wide_in) / 2
  number_of_loci_in_per_sample <- rep(number_of_loci_in, number_of_samples_in)
  number_of_loci_out_per_sample <- unname(sapply(tabs, nrow))
  expect_equal(number_of_loci_in_per_sample, number_of_loci_out_per_sample)

  # verify a couple of datapoints
  expect_identical(wide_in["S1", "vWA"], tabs$S1$Allele1[2])
  expect_identical(wide_in["S2", "vWA(2)"], tabs$S2$Allele2[2])
})

test_that("Wide references to allele tables roundtrip", {

  wide_in <- structure(list(D3S1358 = c("18", "17", "17", "17"),
                            `D3S1358(2)` = c("15","17", "18", "15"),
                            vWA = c("16", "18", "18", "18"),
                            `vWA(2)` = c("16","18", "16", "16"),
                            D16S539 = c("10", "9", "9", "11"),
                            `D16S539(2)` = c("13","11", "13", "10")),
                       class = "data.frame",
                       row.names = c("1", "2", "S1", "S2"))
  # convert
  tabs_out <- .wide_references_to_allele_tables(wide_in)

  # back
  wide_out <- .allele_tables_to_wide_references(tabs_out)

  expect_identical(wide_in, wide_out)
})
