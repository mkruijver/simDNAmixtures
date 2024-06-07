library(simDNAmixtures)

ped_fs <- pedtools::nuclearPed(nch = 2)

freqs <- list(
  locus_a = setNames(c(0.4, 0.6), nm = letters[1:2]),
  locus_b = setNames(c(0.9, 0.1), nm = letters[1:2])
)

linkage_map <- data.frame(chromosome = rep("1",2),
                                 locus = c("locus_a", "locus_b"),
                                 position = c(1.1,1.1))

test_that("Basic example works without linkage", {

  # sample genotypes based on frequencies but no linkage map
  x <- sample_many_pedigree_genotypes(pedigree = ped_fs,
                                      freqs = freqs)

  # check results are sensible
  expect_identical(nrow(x), 4L)
  expect_identical(ncol(x), 4L)
  expect_true(all(as.vector(x) %in% c("a", "b")))
})


test_that("Basic example works with linkage", {

  # sample genotypes based on frequencies and linkage map
  x <- sample_many_pedigree_genotypes(pedigree = ped_fs,
                                      freqs = freqs,
                                      linkage_map = linkage_map)

  # check results are sensible
  expect_identical(nrow(x), 4L)
  expect_identical(ncol(x), 4L)
  expect_true(all(as.vector(x) %in% c("a", "b")))
})


test_that("Extra unrelated persons can be added", {

  # sample genotypes based on frequencies and linkage map
  x <- sample_many_pedigree_genotypes(pedigree = ped_fs,
                                      freqs = freqs,
                                      number_of_extra_unrelateds = 5,
                                      number_of_replicates = 2,
                                      linkage_map = linkage_map)

  x

  # check results are sensible
  expect_identical(nrow(x), 18L)
  expect_identical(ncol(x), 4L)
  expect_true(all(as.vector(x) %in% c("a", "b")))
})

