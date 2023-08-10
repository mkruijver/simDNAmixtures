#' @title Sample log normal model(s) with parameters according to priors
#'
#' @param number_of_contributors Integer
#' @param sampling_parameters List. Needs to contain:
#' \itemize{
#'  \item min_template. Numeric of length one.
#'  \item max_template. Numeric of length one.
#'  \item degradation_shape. Numeric of length one.
#'  \item degradation_scale. Numeric of length one.
#' }
#' @param model_settings List. See \link{log_normal_model}.
#' @details In simulation studies involving many mixed DNA profiles, one often needs to generate various samples with different model parameters. This function samples a log normal model with parameters according to prior distributions. The template parameter for each contributor is sampled uniformly between \code{min_template} and \code{max_template}. The degradation parameter for each contributor is sampled from a gamma distribution with parameters \code{degradation_shape} and \code{degradation_scale}.
#' @return When \code{length(number_of_contributors)==1}, a single \link{log_normal_model} of class \code{pg_model}. Otherwise, a list of these.
#' @examples
#' data(gf)
#'
#' sampling_parameters <- list(min_template = 50., max_template = 1000.,
#'                             degradation_shape = 2.5, degradation_scale = 1e-3)
#'
#' model_no_stutter <- sample_log_normal_model(number_of_contributors = 1,
#'                                             sampling_parameters = sampling_parameters,
#'                                             model_settings = gf$log_normal_settings)
#' @export
sample_log_normal_model <- function(number_of_contributors, sampling_parameters, model_settings){

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

  min_template <- sampling_parameters$min_template
  max_template <- sampling_parameters$max_template
  degradation_shape <- sampling_parameters$degradation_shape
  degradation_scale <- sampling_parameters$degradation_scale

  if ((!is.numeric(min_template)) |
      ((length(min_template) != 1) && (length(min_template) != number_of_contributors))){
    stop("min_template needs to be a numeric of length 1 or number_of_contributors")
  }
  if ((!is.numeric(max_template)) |
      ((length(max_template) != 1) && (length(max_template) != number_of_contributors))){
    stop("max_template needs to be a numeric of length 1 or number_of_contributors")
  }

  if (length(min_template) == 1) min_template <- rep(min_template, number_of_contributors)
  if (length(max_template) == 1) max_template <- rep(max_template, number_of_contributors)

  if (any(min_template <= 0)){
    stop("min_template needs to be positive")
  }
  if (any(max_template <= 0)){
    stop("max_template needs to be positive")
  }
  if (any(min_template > max_template)){
    stop("min_template > max_template")
  }

  validate_log_normal_model_settings(model_settings,
                                     validate_k2 = FALSE,
                                     validate_LSAE = FALSE)

  template <- stats::runif(n = number_of_contributors,
                         min = min_template, max = max_template)

  degradation <- pmin(stats::rgamma(n = number_of_contributors,
                             shape = degradation_shape,
                             scale = degradation_scale),
                      model_settings$degradation_parameter_cap)

  c2_prior <- model_settings$c2_prior
  c2 <- stats::rgamma(n = 1, shape = c2_prior[1],
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
