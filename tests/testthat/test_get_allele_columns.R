library(simDNAmixtures)

test_that("get-allele-columns works for 2 columns", {

  freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc_022024.csv",
                             package = "simDNAmixtures"))


  genotypes <-  sample_contributor_genotypes(contributors = c("U1", "U2"),
                                             freqs)

  allele_columns <- .get_allele_columns(genotypes$U1)

  # verify dimensions
  expect_equal(ncol(allele_columns), 2L)
  expect_equal(nrow(allele_columns), nrow(genotypes$U1))
  }
)
