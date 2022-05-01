#' @title Sample log normal stutter variance parameters according to priors
#'
#' @param log_normal_stutter_variability List of variability parameters. See \link{gf} for an example.
#' @return Named numeric with stutter variance parameter for all stutter types. Names are \code{k2} concatenated with the name of the stuter type. See example.
#' @examples
#'  data(gf)
#'  log_normal_stutter_variability <- gf$log_normal_settings$stutter_variability
#'  k2 <- sample_log_normal_stutter_variance(log_normal_stutter_variability)
#' @export
sample_log_normal_stutter_variance <- function(log_normal_stutter_variability){

  if (!is.list(log_normal_stutter_variability)){
    stop("log_normal_stutter_variability is not a list")
  }
  if (length(names(log_normal_stutter_variability)) != length(log_normal_stutter_variability)){
    stop("log_normal_stutter_variability is not a named list")
  }

  for (stutter_name in names(log_normal_stutter_variability)){
    if (is.null(log_normal_stutter_variability[[stutter_name]]$k2_prior)){
      stop("no k2 variable defined for stutter ", stutter_name)
    }
  }

  k2 <- sapply(log_normal_stutter_variability, function(stutter){
    stats::rgamma(n = 1, shape = stutter$k2_prior[1], scale = stutter$k2_prior[2])})

  names(k2) <- paste0("k2", names(log_normal_stutter_variability))

  k2
}
