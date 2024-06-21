library(simDNAmixtures)

test_that("Sampling a few mixtures and writing to disk",{

  freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc_022024.csv",
                                         package = "simDNAmixtures"))
  gf <- gf_configuration()

  sampling_parameters <- list(min_mu = 50., max_mu = 5e3,
                              min_cv = 0.05, max_cv = 0.35,
                              degradation_shape1 = 10, degradation_shape2 = 1)

  temp_results_dir <- withr::local_tempdir()

  mixtures <- sample_mixtures(n = 2, contributors = c("U1", "U2"), freqs = freqs,
                              sampling_parameters = sampling_parameters,
                              model_settings = gf$gamma_settings_no_stutter,
                              sample_model = sample_gamma_model,
                              results_directory = temp_results_dir,
                              silent = TRUE)
  ## TODO: add verification

  expect_equal(length(mixtures$samples), 2)
})

test_that("Sampling log-normal mixtures with AMEL",{

  freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc_022024.csv",
                                         package = "simDNAmixtures"))
  gf <- gf_configuration()

  sampling_parameters <- list(min_template = 150., max_template = 2000.,
                              degradation_shape = 2.5, degradation_scale = 1e-3)

  # temp_results_dir <- withr::local_tempdir()

  mixtures <- sample_mixtures(n = 2, contributors = c("U1", "U2"), freqs = freqs,
                              sampling_parameters = sampling_parameters,
                              model_settings = gf$log_normal_bwfw_settings,
                              sample_model = sample_log_normal_model,
                              silent = TRUE)


  expect_true("AMEL" %in% mixtures$samples[[1]]$mixture$Locus)
  expect_false(anyNA(mixtures$samples[[1]]$mixture$Height))
  expect_false(anyNA(mixtures$samples[[1]]$mixture$Size))
})
