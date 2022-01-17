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


  gf$stutter_model <- allele_specific_stutter_model(stutter_types = gf$stutters,
                                                    size_regression = gf$size_regression)

  # log-normal stutter variability model
  gf$log_normal_stutter_variability <- list(
    BackStutter = list(k2_prior = c(1.884,7.686),
                       stutter_max = 0.3,
                       inversely_proportional_to_parent = TRUE),
    ForwardStutter = list(k2_prior = c(2.144,4.507),
                       stutter_max = 0.15,
                       inversely_proportional_to_parent = FALSE)
  )

  gf$log_normal_LSAE_variance <- 0.019

  gf
}
