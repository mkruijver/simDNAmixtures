#' Stutters and size regressions for a GlobalFiler 3500 kit
#'
#' A dataset containing default parameters and settings for a GlobalFiler 3500 kit.
#'
#' @format A list of:
#'  \describe{
#'   \item{autosomal_markers}{Names of autosomal markers in the GlobalFiler kit}
#'   \item{repeat_length_by_marker}{Named numeric with STR repeat length by locus name}
#'   \item{size_regression}{See \link{read_size_regression}}
#'   \item{stutters}{List of 4 stutter types, to be used with \link{allele_specific_stutter_model}}
#'   \item{stutter_model}{For convenience, a pre-defined allele_specific_stutter_model}
#'   \item{log_normal_settings}{Settings corresponding to a log normal model with all stutter types}
#'   \item{log_normal_settings_fwbw}{Settings corresponding to a log normal model with backward and forward stutter only}
#'   \item{gamma_settings}{Settings corresponding to a gamma model with all stutter types}
#'   \item{gamma_settings_no_stutter}{Settings for a gamma model without stutter}
#' }
#' @export
gf_configuration <- function(){

  gf <- list()

  # load kits data
  local_env <- environment()
  utils::data(kits, envir = local_env)
  kits <- local_env$kits

  gf$autosomal_markers <- unique(kits$GlobalFiler_Panel_v1$Marker[!kits$GlobalFiler_Panel_v1$Gender.Marker])
  loci_with_AMEL <- unique(kits$GlobalFiler_Panel_v1$Marker[(!kits$GlobalFiler_Panel_v1$Gender.Marker) |
                                                              (kits$GlobalFiler_Panel_v1$Marker == "AMEL")])

  repeat_length_by_marker <- stats::setNames(kits$GlobalFiler_Panel_v1$Repeat[
    match(gf$autosomal_markers, kits$GlobalFiler_Panel_v1$Marker)],
    gf$autosomal_markers)

  gf$repeat_length_by_marker <- repeat_length_by_marker

  # size regression
  filename_size_regression <- system.file("extdata","GlobalFiler_SizeRegression.csv",
                                          package = "simDNAmixtures")
  gf$size_regression <- read_size_regression(filename_size_regression,
                        exceptions = list(AMEL = setNames(c(98.5, 104.5), nm = c("X", "Y"))))

  gf$stutters <- list()

  # back stutter
  filename_bs_exceptions <- system.file("extdata","GlobalFiler_Stutter_Exceptions_3500.csv",
                                        package = "simDNAmixtures")
  bs_exceptions <- read_stutter_exceptions(filename_bs_exceptions)

  filename_bs_regression <- system.file("extdata","GlobalFiler_Stutter_3500.txt",
                                        package = "simDNAmixtures")
  bs_regression <- read_stutter_regression(filename_bs_regression)
  back_stutter <- stutter_type(name = "BackStutter", delta = -1,
                               stutter_regression = bs_regression,
                               stutter_exceptions = bs_exceptions)

  gf$stutters$BackStutter <- back_stutter

  # forward stutter
  filename_fs_regression <- system.file("extdata","GlobalFiler_Forward_Stutter_3500.txt",
                                        package = "simDNAmixtures")
  fs_regression <- read_stutter_regression(filename_fs_regression)
  forward_stutter <- stutter_type(name = "ForwardStutter", delta = 1,
                                  stutter_regression = fs_regression)

  gf$stutters$ForwardStutter <- forward_stutter

  # 2bp back stutter
  filename_2bp_regression <- system.file("extdata","GlobalFiler_2bp_Stutter_3500.txt",package = "simDNAmixtures")
  regression_2bp <- read_stutter_regression(filename_2bp_regression)
  gf$stutters[["2bpBackStutter"]] <- stutter_type(name = "2bpBackStutter",
                                                  delta = c(0, -2),
                                                  repeat_length_by_marker = repeat_length_by_marker,
                                                  applies_to_all_loci = FALSE,
                                                  applies_to_loci = c("SE33", "D1S1656"),
                                                  stutter_regression = regression_2bp)

  # double back stutter
  filename_double_bs_regression <- system.file("extdata","GlobalFiler_Double_Back_Stutter_3500.txt",
                                               package = "simDNAmixtures")
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

  detection_threshold <- c(D3S1358 = 75, vWA = 75, D16S539 = 75, CSF1PO = 75, TPOX = 75, AMEL = 75,
                           D8S1179 = 100, D21S11 = 100, D18S51 = 100, D2S441 = 60, D19S433 = 60,
                           TH01 = 60, FGA = 60, D22S1045 = 80, D5S818 = 80, D13S317 = 80,
                           D7S820 = 80, SE33 = 80, D10S1248 = 100, D1S1656 = 100, D12S391 = 100,
                           D2S1338 = 100)

  gf$log_normal_settings <- list(
    locus_names = loci_with_AMEL,
    degradation_parameter_cap = 0.01,
    c2_prior = c(8.45,1.746),
    LSAE_variance_prior = 0.019,
    detection_threshold = detection_threshold,
    size_regression = gf$size_regression,
    stutter_model = gf$stutter_model,
    stutter_variability = log_normal_stutter_variability
  )

  stutter_model_bwfw <- allele_specific_stutter_model(stutter_types = gf$stutters[c(
    "BackStutter", "ForwardStutter")],
    size_regression = gf$size_regression)

  detection_threshold_75 <- c(D3S1358 = 75, vWA = 75, D16S539 = 75, CSF1PO = 75, TPOX = 75, AMEL = 75,
                              D8S1179 = 75, D21S11 = 75, D18S51 = 75, D2S441 = 75, D19S433 = 75,
                              TH01 = 75, FGA = 75, D22S1045 = 75, D5S818 = 75, D13S317 = 75,
                              D7S820 = 75, SE33 = 75, D10S1248 = 75, D1S1656 = 75, D12S391 = 75,
                              D2S1338 = 75)

  log_normal_stutter_variability_bwfw <- list(
    BackStutter = list(k2_prior = c(3.499, 4.803),
                       inversely_proportional_to_parent = TRUE,
                       max_stutter_ratio = 0.3),
    ForwardStutter = list(k2_prior = c(4.865, 3.101),
                          inversely_proportional_to_parent = FALSE,
                          max_stutter_ratio = 0.15))


  gf$log_normal_bwfw_settings <- list(
    locus_names = loci_with_AMEL,
    degradation_parameter_cap = 0.01,
    c2_prior = c(4.865, 3.101),
    LSAE_variance_prior = 0.0217,
    detection_threshold = detection_threshold_75,
    size_regression = gf$size_regression,
    stutter_model = stutter_model_bwfw,
    stutter_variability = log_normal_stutter_variability_bwfw
  )

  gf$gamma_settings <- list(
    locus_names = loci_with_AMEL,
    detection_threshold = detection_threshold,
    LSAE_variance_prior = 0,
    size_regression = gf$size_regression,
    stutter_model = gf$stutter_model
  )

  gf$gamma_settings_no_stutter <- list(
    locus_names = loci_with_AMEL,
    detection_threshold = detection_threshold,
    LSAE_variance_prior = 0,
    size_regression = gf$size_regression
  )

  gf
}
