#' @title Sample drop model(s) with parameters according to priors
#'
#' @param number_of_contributors Integer
#' @param sampling_parameters List. Needs to contain:
#' \itemize{
#'  \item min_dropout_probability. Numeric of length one.
#'  \item max_dropout_probability Numeric of length one.
#' }
#' @param drop_in_rate Numeric vector of length one. Expected number of drop-ins per locus. Default is 0.
#' @param model_settings List. See \link{drop_model}.
#' @details In simulation studies involving many mixed DNA profiles, one often needs to generate various samples with different model parameters. This function samples a drop model with parameters according to prior distributions. The dropout probability for each contributor is sampled uniformly between \code{min_dropout_probability.} and \code{max_dropout_probability}.
#' @return When \code{length(number_of_contributors)==1}, a single \link{drop_model} of class \code{pg_model}. Otherwise, a list of these.
#' @seealso [sample_mixtures_fixed_parameters] to directly supply parameters of choice for more control
#' @examples
#' gf <- gf_configuration()
#'
#' sampling_parameters <- list(min_dropout_probability. = 0., max_dropout_probability. = 0.5)
#'
#' model <- sample_drop_model(number_of_contributors = 1,
#'                            sampling_parameters = sampling_parameters,
#'                            model_settings = list(locus_names = gf$autosomal_markers,
#'                            size_regression = gf$size_regression))
#' @export
sample_drop_model <- function(number_of_contributors, sampling_parameters,
                              drop_in_rate = 0.,
                              model_settings){

  if (length(number_of_contributors) > 1){
    this_call <- match.call()

    return(lapply(number_of_contributors,
                  function(n){
                    this_call$number_of_contributors <- n
                    eval(this_call)
                  }))
  }


  d <- stats::runif(n = number_of_contributors,
                    min = sampling_parameters$min_dropout_probability,
                    max = sampling_parameters$max_dropout_probability)

  model <- drop_model(dropout_probabilities = d,
                      drop_in_rate = drop_in_rate,
                      model_settings = model_settings)

  model
}
