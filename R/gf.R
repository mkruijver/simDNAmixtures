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
"gf"
