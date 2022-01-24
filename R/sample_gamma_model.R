#' @title Sample gamma model(s) with parameters according to priors
#'
#' @param number_of_contributors Integer
#' @param sampling_parameters List. Needs to contain:
#' \itemize{
#'  \item min_mu. Numeric of length one.
#'  \item max_mu. Numeric of length one.
#'  \item min_cv. Numeric of length one.
#'  \item max_cv. Numeric of length one.
#'  \item degradation_shape1. Numeric of length one.
#'  \item degradation_shape2. Numeric of length one.
#' }
#' @param model_settings List. See \link{gamma_model}.
#' @examples
#' gf <- get_GlobalFiler_3500_data()
#'
#' sampling_parameters <- list(min_mu = 50., max_mu = 5e3,
#'                            min_cv = 0.05, max_cv = 0.35,
#'                            degradation_shape1 = 10, degradation_shape2 = 1)
#'
#' model_no_stutter <- sample_gamma_model(number_of_contributors = 2,
#'                                       sampling_parameters = sampling_parameters,
#'                                       model_settings = gf$gamma_settings_no_stutter)
#'
#' model_no_stutter$parameters
#'
#' @export
sample_gamma_model <- function(number_of_contributors, sampling_parameters, model_settings){

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

  if (as.character(number_of_contributors) != as.character(as.integer(number_of_contributors))){
    stop("number_of_contributors needs to be integer valued")
  }

  number_of_contributors <- as.integer(number_of_contributors)

  min_mu <- sampling_parameters$min_mu
  max_mu <- sampling_parameters$max_mu

  min_cv <- sampling_parameters$min_cv
  max_cv <- sampling_parameters$max_cv

  degradation_shape1 <- sampling_parameters$degradation_shape1
  degradation_shape2 <- sampling_parameters$degradation_shape2

  if ((!is.numeric(min_mu)) | (length(min_mu) != 1)){
    stop("min_mu needs to be a numeric of length 1")
  }
  if ((!is.numeric(max_mu)) | (length(max_mu) != 1)){
    stop("max_mu needs to be a numeric of length 1")
  }
  if (min_mu <= 0){
    stop("min_mu needs to be positive")
  }
  if (max_mu <= 0){
    stop("max_mu needs to be positive")
  }
  if (min_mu > max_mu){
    stop("min_mu > max_mu")
  }

  if ((!is.numeric(min_cv)) | (length(min_cv) != 1)){
    stop("min_cv needs to be a numeric of length 1")
  }
  if ((!is.numeric(max_cv)) | (length(max_cv) != 1)){
    stop("max_cv needs to be a numeric of length 1")
  }
  if (min_cv <= 0){
    stop("min_cv needs to be positive")
  }
  if (max_cv <= 0){
    stop("max_cv needs to be positive")
  }
  if (min_cv > max_cv){
    stop("min_cv > max_cv")
  }

  mu <- runif(n = 1, min = min_mu, max = max_mu)

  cv <- runif(n = 1, min = min_cv, max = max_cv)

  mixture_proportions_unnormalised <- sort(
              runif(n = number_of_contributors,
                    min = 0, max = 1),
              decreasing = TRUE)

  mixture_proportions <- mixture_proportions_unnormalised /
    sum(mixture_proportions_unnormalised)

  degradation <- rbeta(n = number_of_contributors,
                       shape1 = degradation_shape1, shape2 = degradation_shape2)

  LSAE <- sample_LSAE(model_settings$LSAE_variance_prior,
                      model_settings$locus_names)

  model <- gamma_model(mixture_proportions = mixture_proportions, mu = mu,
                       cv = cv, degradation_beta = degradation, LSAE = LSAE,
                       model_settings = model_settings)

  model
}
