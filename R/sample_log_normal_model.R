#' @title Sample log normal model(s) with parameters according to priors
#'
#' @param number_of_contributors Integer
#' @param min_template Numeric of length one.
#' @param max_template Numeric of length one.
#' @param degradation_shape Numeric of length one.
#' @param degradation_scale Numeric of length one.
#' @param model_settings List. Possible parameters: \itemize{
#'  \item locus_names. Character vector.
#'  \item degradation_parameter_cap. Numeric.
#'  \item c2_prior. Numeric of length two with shape and scale.
#'  \item LSAE_variance_prior. Numeric of length one.
#'  \item detection_threshold. Numeric vector (named) with Detection Thresholds. Defaults to 50 for each locus.
#'  \item size_regression. Function, see \link{read_size_regression}.
#'  \item stutter_model. Optionally a stutter_model object that gives expected stutter heights. See \link{global_stutter_model}.
#'  \item stutter_variability. Optionally peak height variability parameters for stutters. Required when stutter_model is supplied.
#' }
#' @examples
#' gf <- get_GlobalFiler_3500_data()
#' model_no_stutter <- sample_log_normal_model(number_of_contributors = 1,
#'                                             min_template = 50., max_template = 1000.,
#'                                             degradation_shape = 2.5, degradation_scale = 1e-3,
#'                                             model_settings = gf$log_normal_settings)
#' @export
sample_log_normal_model <- function(number_of_contributors,
                                    min_template, max_template,
                                    degradation_shape, degradation_scale,
                                    model_settings){

  if (length(number_of_contributors) > 1){
    this_call <- match.call()

    return(lapply(number_of_contributors,
                  function(n){
                    this_call$number_of_contributors <- n
                    eval(this_call)
                  }))
  }

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

  validate_log_normal_model_settings(model_settings,
                                     validate_k2 = FALSE,
                                     validate_LSAE = FALSE)

  template <- sort(runif(n = number_of_contributors,
                         min = min_template, max = max_template),
                   decreasing = TRUE)

  degradation <- pmin(rgamma(n = number_of_contributors,
                             shape = degradation_shape,
                             scale = degradation_scale),
                      model_settings$degradation_parameter_cap)

  c2_prior <- model_settings$c2_prior
  c2 <- rgamma(n = 1, shape = c2_prior[1],
               scale = c2_prior[2])

  LSAE <- sample_LSAE(model_settings$LSAE_variance_prior,
                      model_settings$locus_names)

  if (is.null(model_settings$stutter_model)){
    model <- log_normal_model(template = template, degradation = degradation,
                              LSAE = LSAE, c2 = c2,
                              model_settings = model_settings)
  }
  else{
    k2 <- sample_log_normal_stutter_variance(model_settings$stutter_variability)

    model <- log_normal_model(template = template, degradation = degradation,
                              LSAE = LSAE, c2 = c2,
                              k2 = k2,
                              model_settings = model_settings)
  }

  model
}
