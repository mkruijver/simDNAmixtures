
test_that("Wide YSTR references (df) convert to allele tables", {

  x <- structure(list(DYS576 = c("18", "17"), DYS389I = c("13", "13"),
                      DYS635 = c("23", "23"), DYS389II = c("29", "29"),
                      DYS627 = c("22", "21"), DYS460 = c("11", "10"),
                      DYS458 = c("17", "17"), DYS19 = c("14", "14"),
                      YGATAH4 = c("11", "12"),DYS448 = c("19", "19"),
                      DYS391 = c("11", "11"), DYS456 = c("16", "15"),
                      DYS390 = c("24", "23"), DYS438 = c("12", "12"),
                      DYS392 = c("13", "13"), DYS518 = c("41", "40"),
                      DYS570 = c("18", "18"), DYS437 = c("15", "15"),
                      DYS385 = c("12,13", "11,14"), DYS449 = c("29", "28"),
                      DYS393 = c("13", "13"), DYS439 = c("11", "12"),
                      DYS481 = c("22", "22"), DYF387S1 = c("34,36", "35,36"),
                      DYS533 = c("12", "11")),
                 row.names = c("K1", "K2"), class = "data.frame")

  # convert
  tabs <- .wide_YSTR_references_to_allele_tables(x)

  # verify dimensions
  number_of_samples_in <- nrow(x)
  number_of_samples_out <- length(tabs)
  expect_identical(number_of_samples_in, number_of_samples_out)

  number_of_loci_in <- ncol(x)
  number_of_loci_in_per_sample <- rep(number_of_loci_in, number_of_samples_in)
  number_of_loci_out_per_sample <- unname(sapply(tabs, nrow))
  expect_equal(number_of_loci_in_per_sample, number_of_loci_out_per_sample)

  # verify a couple of datapoints
  expect_identical(x["K1", "DYS627"], tabs$K1$Allele1[tabs$K2$Locus=="DYS627"])
  expect_identical(x["K2", "DYS456"], tabs$K2$Allele1[tabs$K2$Locus=="DYS456"])

  expect_equal(x["K2", "DYF387S1"],
                paste0(unname(tabs$K2[tabs$K2$Locus=="DYF387S1",][3:4]), collapse = ","))
})
