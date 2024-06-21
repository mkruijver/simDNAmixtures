test_that("Sampling founders without pedigree", {
  freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc_022024.csv",
                                         package = "simDNAmixtures"))

  unrelated_names <- c("K1", "K2")
  x <- .sample_many_founders(freqs = freqs, unrelated_names = unrelated_names,
                             loci = c("AMEL", "vWA"))

  expect_equal(rownames(x), unrelated_names)
})

test_that("Sampling founders works with AMEL", {

  freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc_022024.csv",
                                         package = "simDNAmixtures"))

  # sample ped 1 (father, mother, son1, son2)
  ped1 <- pedtools::nuclearPed(nch = 2, sex = c(1, 1))
  x1 <- .sample_many_founders(ped = ped1, freqs = freqs, loci = c("AMEL", "vWA"))
  expect_identical(nrow(x1), 4L)
  expect_identical(ncol(x1), 4L)

  expect_equal(as.vector(x1["3", c("AMEL", "AMEL(2)")]), c("X", "Y"))
  expect_equal(as.vector(x1["4", c("AMEL", "AMEL(2)")]), c("X", "Y"))

  # sample ped 2 (father, mother, son1, daughter 1)
  ped2 <- pedtools::nuclearPed(nch = 2, sex = c(1, 2))
  x2 <- .sample_many_founders(ped = ped2, freqs = freqs, loci = c("AMEL", "vWA"))
  expect_identical(nrow(x2), 4L)
  expect_identical(ncol(x2), 4L)

  expect_equal(as.vector(x2["3", c("AMEL", "AMEL(2)")]), c("X", "Y"))
  expect_equal(as.vector(x2["4", c("AMEL", "AMEL(2)")]), c("X", "X"))

  # sample ped 3 (father, mother, daughter 1, 2 extra unrelateds)
  ped3 <- pedtools::nuclearPed(nch = 1, sex = c(2))
  x3 <- .sample_many_founders(ped = ped3, freqs = freqs, loci = c("AMEL", "vWA"),
                              unrelated_names = c("U1", "U2"))
  expect_identical(nrow(x3), 5L)
  expect_identical(ncol(x3), 4L)

  expect_equal(as.vector(x3["3", c("AMEL", "AMEL(2)")]), c("X", "X"))

  # reps
  x3_reps <- .sample_many_founders(ped = ped3, freqs = freqs, loci = c("AMEL", "vWA"),
                                   unrelated_names = c("U1", "U2"),
                                   number_of_replicates = 10)

  expect_identical(nrow(x3_reps), 5L * 10L)
  expect_identical(ncol(x3), 4L)
})


test_that("Sampling founders without pedigree with AMEL and other loci", {
  freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc_022024.csv",
                                         package = "simDNAmixtures"))

  unrelated_names <- c("K1", "K2")
  x <- .sample_many_founders(freqs = freqs, unrelated_names = unrelated_names,
                             loci = c("vWA", "TH01", "AMEL"))

  expect_false(anyNA(x))
})
