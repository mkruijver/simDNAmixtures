test_that("Sampling contributor genotypes without pedigree", {
  freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc_022024.csv",
                                         package = "simDNAmixtures"))

  unrelated_names <- c("K1", "K2")
  x <- .sample_many_founders(freqs = freqs, unrelated_names = unrelated_names,
                             loci = c("AMEL", "vWA"))

  expect_equal(rownames(x), unrelated_names)
})
