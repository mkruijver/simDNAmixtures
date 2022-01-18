#' @title Sample log normal model with parameters according to priors
#'
#' @param number_of_contributors Integer
#' @param min_template Numeric of length one.
#' @param max_template Numeric of length one.
#' @param degradation_shape Numeric of length one.
#' @param degradation_scale Numeric of length one.
#' @param locus_names Character vector.
#' @param degradation_parameter_cap Numeric. Defaults to 0.01.
#' @param c2_prior Numeric of length two with shape and scale.
#' @param LSAE_variance_prior Numeric of length one.
#' @param detection_threshold Numeric vector (named) with Detection Thresholds. Defaults to 50 for each locus.
#' @param size_regression Function, see \link{read_size_regression}.
#' @param stutter_model Optionally a stutter_model object that gives expected stutter heights. See \link{global_stutter_model}.
#' @param stutter_variability Optionally peak height variability parameters for stutters. Required when stutter_model is supplied.
#' @examples
#' gf <- get_GlobalFiler_3500_data()
#' model_no_stutter <- sample_log_normal_model(number_of_contributors = 1,
#'                                             min_template = 50., max_template = 1000.,
#'                                             degradation_shape = 2.5, degradation_scale = 1e-3,
#'                                             locus_names = gf$autosomal_markers,
#'                                             c2_prior = gf$log_normal_c2_prior,
#'                                             LSAE_variance_prior = gf$log_normal_LSAE_variance,
#'                                             size_regression = gf$size_regression)
#' @export
sample_log_normal_model <- function(number_of_contributors,
                                    min_template, max_template,
                                    degradation_shape, degradation_scale,
                                    locus_names,
                                    degradation_parameter_cap = 0.01,
                                    c2_prior,
                                    LSAE_variance_prior,
                                    detection_threshold = setNames(rep(50., length(locus_names)), locus_names),
                                    size_regression,
                                    stutter_model,
                                    stutter_variability){

  if (length(number_of_contributors) != 1){
    stop("number_of_contributors needs to have length 1")
  }

  if (!(is.numeric(number_of_contributors) | is.integer(number_of_contributors))){
    stop("number_of_contributors needs to be integer valued")
  }

  if (as.character(number_of_contributors)!=as.character(as.integer(number_of_contributors))){
    stop("number_of_contributors needs to be integer valued")
  }

  number_of_contributors <- as.integer(number_of_contributors)

  if ((!is.numeric(min_template)) | (length(min_template) != 1)){
    stop("min_template needs to be a numeric of length 1")
  }
  if ((!is.numeric(max_template)) | (length(max_template) != 1)){
    stop("max_template needs to be a numeric of length 1")
  }
  if (min_template <= 0){
    stop("min_template needs to be positive")
  }
  if (max_template <= 0){
    stop("max_template needs to be positive")
  }
  if (min_template > max_template){
    stop("min_template > max_template")
  }
  if ((!is.numeric(degradation_parameter_cap)) | (length(degradation_parameter_cap) != 1)){
    stop("degradation_parameter_cap needs to be a numeric of length 1")
  }

  if (!missing(stutter_model)){
    if (missing(stutter_variability)){
      stop("stutter_variability needs to be supplied when stutter_model is supplied")
    }

    for (stutter_name in names(stutter_variability)){
      if (!is.numeric(stutter_variability[[stutter_name]]$max_stutter_ratio)){
        stop("stutter_variability$", stutter_name,
             "$max_stutter_ratio is not numeric")
      }
      if (length(stutter_variability[[stutter_name]]$max_stutter_ratio) != 1){
        stop("stutter_variability$", stutter_name,
             "$max_stutter_ratio needs to have length 1")
      }
      if (stutter_variability[[stutter_name]]$max_stutter_ratio < 0){
        stop("stutter_variability$", stutter_name,
             "$max_stutter_ratio needs to be non-negative")
      }
    }
  }

  template <- sort(runif(n = number_of_contributors,
                         min = min_template, max = max_template),
                   decreasing = TRUE)

  degradation <- pmin(rgamma(n = number_of_contributors,
                             shape = degradation_shape,
                             scale = degradation_scale),
                      degradation_parameter_cap)

  c2 <- rgamma(n = 1, shape = c2_prior[1],
               scale = c2_prior[2])

  LSAE <- sample_LSAE(LSAE_variance_prior, locus_names)

  if (missing(stutter_model)){
    model <- log_normal_model(template = template, degradation = degradation,
                              locus_names = locus_names, LSAE = LSAE, c2 = c2,
                              size_regression = size_regression)
  }
  else{
    k2 <- sample_log_normal_stutter_variance(stutter_variability)

    model <- log_normal_model(template = template, degradation = degradation,
                              locus_names = locus_names, LSAE = LSAE, c2 = c2,
                              k2 = k2, stutter_model = stutter_model,
                              stutter_variability = stutter_variability,
                              size_regression = size_regression)
  }

  model
}
