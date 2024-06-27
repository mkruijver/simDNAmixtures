#' K1 <- data.frame(
#'   `Sample Name` = "K1",
#'   Locus = c("D3S1358", "vWA", "D16S539", "CSF1PO", "TPOX", "D8S1179",
#'             "D21S11", "D18S51", "D2S441", "D19S433", "TH01", "FGA",
#'             "D22S1045", "D5S818", "D13S317", "D7S820", "SE33",
#'             "D10S1248", "D1S1656", "D12S391", "D2S1338"),
#'   Allele1 = c("14", "15", "11", "11", "8",
#'               "11", "28", "14", "11", "14", "9.3",
#'               "21", "15", "9", "9", "10", "14",
#'               "14", "12", "21", "20"),
#'   Allele2 = c("16", "18", "13", "11", "11",
#'               "12", "30", "16", "12", "14", "9.3",
#'               "24", "16", "11", "12", "10", "18",
#'               "14", "15", "22", "24"), check.names = FALSE
#' )
#'
#' K2 <- data.frame(
#'   `Sample Name` = "K2",
#'   Locus = c("D3S1358", "vWA", "D16S539", "CSF1PO", "TPOX", "D8S1179",
#'             "D21S11", "D18S51", "D2S441", "D19S433", "TH01", "FGA",
#'             "D22S1045", "D5S818", "D13S317", "D7S820", "SE33",
#'             "D10S1248", "D1S1656", "D12S391", "D2S1338"),
#'   Allele1 = c("16", "16", "10", "13", "8",
#'               "11", "27", "17", "10", "13", "6",
#'               "23", "16", "11", "11", "9", "22.2",
#'               "14", "12", "22", "21"),
#'   Allele2 = c("16", "17", "12", "14", "8",
#'               "13", "30", "18", "11", "14", "6",
#'               "25", "17", "11", "11", "12", "28.2",
#'               "15", "14", "22", "23"), check.names = FALSE
#' )
#'
#' # first sample three replicates of a low-level profile of K1 only
#' gf <- gf_configuration()
#'
#' sampling_parameters <- list(min_template = 75., max_template = 75,
#'                             degradation_shape = 2.5, degradation_scale = 1e-3)
#'
#' single_source_results <- sample_mixtures_from_genotypes(n = 1,
#'                 genotypes = list(K1), sampling_parameters = sampling_parameters,
#'                 number_of_replicates = 3, sample_model = sample_log_normal_model,
#'                 model_settings = gf$log_normal_bwfw_settings)
#'


test_that("Drop model sampling (single source, known reference)", {

  K1 <- data.frame(
    `Sample Name` = "K1",
    Locus = c("D3S1358", "vWA", "D16S539", "CSF1PO", "TPOX", "D8S1179",
              "D21S11", "D18S51", "D2S441", "D19S433", "TH01", "FGA",
              "D22S1045", "D5S818", "D13S317", "D7S820", "SE33",
              "D10S1248", "D1S1656", "D12S391", "D2S1338"),
    Allele1 = c("14", "15", "11", "11", "8",
                "11", "28", "14", "11", "14", "9.3",
                "21", "15", "9", "9", "10", "14",
                "14", "12", "21", "20"),
    Allele2 = c("16", "18", "13", "11", "11",
                "12", "30", "16", "12", "14", "9.3",
                "24", "16", "11", "12", "10", "18",
                "14", "15", "22", "24"), check.names = FALSE
  )

  gf <- gf_configuration()


  freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc_022024.csv",
                                         package = "simDNAmixtures"))

  model_settings <-list(locus_names = gf$autosomal_markers,
                        size_regression = gf$size_regression)

  model <- drop_model(dropout_probabilities = 0., freqs = freqs,
                      model_settings = model_settings)

  x <- sample_mixture_from_genotypes(genotypes = list(K1), model)

  # verify that nothing dropped out
  sample_marker_allele <- paste0(x$Marker, "_", x$Allele)

  ref_marker_allele <- c(paste0(K1$Locus, "_", K1$Allele1),
                         paste0(K1$Locus, "_", K1$Allele2))

  expect_setequal(sample_marker_allele, ref_marker_allele)
})

test_that("Drop model sampling (2P, known references)", {

  K1 <- data.frame(
    `Sample Name` = "K1",
    Locus = c("D3S1358", "vWA", "D16S539", "CSF1PO", "TPOX", "D8S1179",
              "D21S11", "D18S51", "D2S441", "D19S433", "TH01", "FGA",
              "D22S1045", "D5S818", "D13S317", "D7S820", "SE33",
              "D10S1248", "D1S1656", "D12S391", "D2S1338"),
    Allele1 = c("14", "15", "11", "11", "8",
                "11", "28", "14", "11", "14", "9.3",
                "21", "15", "9", "9", "10", "14",
                "14", "12", "21", "20"),
    Allele2 = c("16", "18", "13", "11", "11",
                "12", "30", "16", "12", "14", "9.3",
                "24", "16", "11", "12", "10", "18",
                "14", "15", "22", "24"), check.names = FALSE
  )

  K2 <- data.frame(
    `Sample Name` = "K2",
    Locus = c("D3S1358", "vWA", "D16S539", "CSF1PO", "TPOX", "D8S1179",
              "D21S11", "D18S51", "D2S441", "D19S433", "TH01", "FGA",
              "D22S1045", "D5S818", "D13S317", "D7S820", "SE33",
              "D10S1248", "D1S1656", "D12S391", "D2S1338"),
    Allele1 = c("16", "16", "10", "13", "8",
                "11", "27", "17", "10", "13", "6",
                "23", "16", "11", "11", "9", "22.2",
                "14", "12", "22", "21"),
    Allele2 = c("16", "17", "12", "14", "8",
                "13", "30", "18", "11", "14", "6",
                "25", "17", "11", "11", "12", "28.2",
                "15", "14", "22", "23"), check.names = FALSE
  )

  gf <- gf_configuration()

  freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc_022024.csv",
                                         package = "simDNAmixtures"))

  model_settings <-list(locus_names = gf$autosomal_markers,
                        size_regression = gf$size_regression)

  model <- drop_model(dropout_probabilities = c(0., 0.), freqs = freqs,
                      model_settings = model_settings)

  x <- sample_mixture_from_genotypes(genotypes = list(K1, K2), model)

  # verify that nothing dropped out
  sample_marker_allele <- paste0(x$Marker, "_", x$Allele)

  ref_marker_allele <- c(paste0(K1$Locus, "_", K1$Allele1),
                         paste0(K1$Locus, "_", K1$Allele2),
                         paste0(K2$Locus, "_", K2$Allele1),
                         paste0(K2$Locus, "_", K2$Allele2))

  expect_setequal(sample_marker_allele, ref_marker_allele)
})
