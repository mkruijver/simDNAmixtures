#' @title Sample log normal stutter variance parameters according to priors
#'
#' @param log_normal_stutter_variability List of variability parameters. See \link{get_GlobalFiler_3500_data} for an example.
#' @examples
#'  gf <- get_GlobalFiler_3500_data()
#'  log_normal_stutter_variability <- gf$log_normal_stutter_variability
#'  k2 <- sample_log_normal_stutter_variance(log_normal_stutter_variability)
#' @export
sample_log_normal_stutter_variance <- function(log_normal_stutter_variability){

  if (!is.list(log_normal_stutter_variability)){
    stop("log_normal_stutter_variability is not a list")
  }
  if (is.null(names(log_normal_stutter_variability))){
    stop("log_normal_stutter_variability is not a named list")
  }

  for (stutter_name in names(log_normal_stutter_variability)){
    if (is.null(log_normal_stutter_variability[[stutter_name]]$k2_prior)){
      stop("no k2 variable defined for stutter ", stutter_name)
    }
  }

  k2 <- sapply(log_normal_stutter_variability, function(stutter){
    rgamma(n = 1, shape = stutter$k2_prior[1], scale = stutter$k2_prior[2])})

  names(k2) <- paste0("k2", names(log_normal_stutter_variability))

  k2
}