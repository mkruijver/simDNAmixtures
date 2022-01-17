library(SimMixDNA)

test_that("Log-Normal sampling (single source, no stutter)", {

  gf <- get_GlobalFiler_3500_data()
  freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc.csv",package = "SimMixDNA"))

  model <- log_normal_model(template = 1e3,
                            c2 = 15, locus_names = gf$autosomal_markers,
                            size_regression = gf$size_regression)

  g <- sample_genotype(freqs = freqs, loci = gf$autosomal_markers)

  s <- sample_mixture_profile(list(g), model, sample_name = "mix1")

  expect_equal(s$SampleName, rep("mix1", nrow(s)))
})

test_that("Log-Normal sampling (single source, back and forward stutter)", {

  gf <- get_GlobalFiler_3500_data()
  freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc.csv",package = "SimMixDNA"))

  stutter_types <- list(BackStutter = gf$stutters$BackStutter,
                        ForwardStutter = gf$stutters$ForwardStutter)

  stutter_model <- allele_specific_stutter_model(stutter_types = stutter_types,
                                size_regression = gf$size_regression)

  k2 <- sample_log_normal_stutter_variance(gf$log_normal_stutter_variability)

  model <- log_normal_model(template = 1e3,
                            c2 = 15, k2 = k2,
                            locus_names = gf$autosomal_markers,
                            size_regression = gf$size_regression,
                            stutter_model = stutter_model,
                            stutter_variability = gf$log_normal_stutter_variability)

  g <- sample_genotype(freqs = freqs, loci = gf$autosomal_markers)

  s <- sample_mixture_profile(list(g), model, sample_name = "mix1")

  expect_equal(s$SampleName, rep("mix1", nrow(s)))
  expect_true("ExpectedBackStutter" %in% names(s))
  expect_true("ExpectedForwardStutter" %in% names(s))
  expect_true("ExpectedStutter" %in% names(s))

  expect_equal(s$ExpectedStutter, s$ExpectedBackStutter + s$ExpectedForwardStutter)
})
