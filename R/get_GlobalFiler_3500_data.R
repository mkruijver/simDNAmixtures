#' Convenience function that loads size regression and stutter files
#'
#' @return Object of class list
#' @examples
#' gf <- get_GlobalFiler_3500_data()
#' @export
get_GlobalFiler_3500_data <- function() {

  gf <- list()

  kit_data <- get_kit_data()
  gf$autosomal_markers <- unique(kit_data$GlobalFiler_Panel_v1$Marker[!kit_data$GlobalFiler_Panel_v1$Gender.Marker])

  repeat_length_by_marker <- setNames(kit_data$GlobalFiler_Panel_v1$Repeat[
    match(gf$autosomal_markers, kit_data$GlobalFiler_Panel_v1$Marker)],
    gf$autosomal_markers)

  gf$repeat_length_by_marker <- repeat_length_by_marker

  # size regression
  filename_size_regression <- system.file("extdata","GlobalFiler_SizeRegression.csv",package = "SimMixDNA")
  gf$size_regression <- read_size_regression(filename_size_regression)

  # stutters
  gf$stutters <- list()

  # back stutter
  filename_bs_exceptions <- system.file("extdata","GlobalFiler_Stutter_Exceptions_3500.csv",package = "SimMixDNA")
  bs_exceptions <- read_stutter_exceptions(filename_bs_exceptions)

  filename_bs_regression <- system.file("extdata","GlobalFiler_Stutter_3500.txt",package = "SimMixDNA")
  bs_regression <- read_stutter_regression(filename_bs_regression)
  back_stutter <- stutter_type(name = "BackStutter", delta = -1,
                              stutter_regression = bs_regression,
                              stutter_exceptions = bs_exceptions)

  gf$stutters$BackStutter <- back_stutter

  # forward stutter
  filename_fs_regression <- system.file("extdata","GlobalFiler_Forward_Stutter_3500.txt",package = "SimMixDNA")
  fs_regression <- read_stutter_regression(filename_fs_regression)
  forward_stutter <- stutter_type(name = "ForwardStutter", delta = 1,
                              stutter_regression = fs_regression)

  gf$stutters$ForwardStutter <- forward_stutter

  # 2bp back stutter
  filename_2bp_regression <- system.file("extdata","GlobalFiler_2bp_Stutter_3500.txt",package = "SimMixDNA")
  regression_2bp <- read_stutter_regression(filename_2bp_regression)
  gf$stutters[["2bpBackStutter"]] <- stutter_type(name = "2bpBackStutter",
                                              delta = c(0, -2),
                                              repeat_length_by_marker = repeat_length_by_marker,
                                              applies_to_all_loci = FALSE,
                                              applies_to_loci = c("SE33", "D1S1656"),
                                              stutter_regression = regression_2bp)

  # double back stutter
  filename_double_bs_regression <- system.file("extdata","GlobalFiler_Double_Back_Stutter_3500.txt",
                                               package = "SimMixDNA")
  double_bs_regression <- read_stutter_regression(filename_double_bs_regression)
  gf$stutters$DoubleBackStutter <- stutter_type(name = "DoubleBackStutter",
                                                delta = c(-2),
                                                applies_to_all_loci = TRUE,
                                                stutter_regression = double_bs_regression)

  gf$stutter_model <- allele_specific_stutter_model(stutter_types = gf$stutters,
                                                    size_regression = gf$size_regression)

  # log-normal stutter variability model
  log_normal_stutter_variability <- list(
    BackStutter = list(k2_prior = c(1.884, 7.686),
                       inversely_proportional_to_parent = TRUE,
                       max_stutter_ratio = 0.3),
    ForwardStutter = list(k2_prior = c(2.144, 4.507),
                          inversely_proportional_to_parent = FALSE,
                          max_stutter_ratio = 0.15),
    "2bpBackStutter" = list(k2_prior = c(2.189, 1.431),
                            inversely_proportional_to_parent = FALSE,
                            max_stutter_ratio = 0.1),
    DoubleBackStutter = list(k2_prior = c(3.429, 2.032),
                            inversely_proportional_to_parent = FALSE,
                            max_stutter_ratio = 0.05)
  )

  detection_threshold <- c(D3S1358 = 50, vWA = 50, D16S539 = 50, CSF1PO = 50, TPOX = 50,
                           D8S1179 = 50, D21S11 = 50, D18S51 = 50, D2S441 = 50, D19S433 = 50,
                           TH01 = 50, FGA = 50, D22S1045 = 50, D5S818 = 50, D13S317 = 50,
                           D7S820 = 50, SE33 = 50, D10S1248 = 50, D1S1656 = 50, D12S391 = 50,
                           D2S1338 = 50)

  gf$log_normal_settings <- list(
    locus_names = gf$autosomal_markers,
    degradation_parameter_cap = 0.01,
    c2_prior = c(8.45,1.746),
    LSAE_variance_prior = 0.019,
    detection_threshold = detection_threshold,
    size_regression = gf$size_regression,
    stutter_model = gf$stutter_model,
    stutter_variability = log_normal_stutter_variability
  )

  gf$gamma_settings <- list(
    locus_names = gf$autosomal_markers,
    detection_threshold = detection_threshold,
    LSAE_variance_prior = 0,
    size_regression = gf$size_regression,
    stutter_model = gf$stutter_model
  )

  gf$gamma_settings_no_stutter <- list(
    locus_names = gf$autosomal_markers,
    detection_threshold = detection_threshold,
    LSAE_variance_prior = 0,
    size_regression = gf$size_regression
  )

  gf
}
