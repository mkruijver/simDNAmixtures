#' @title Defines a log normal model for peak height variability
#'
#' @param template Numeric vector
#' @param degradation Numeric Vector of same length as template. Degradation parameters for each contributor.
#' @param c2 Numeric. Allele variance parameter.
#' @param size_regression. Function, see \link{read_size_regression}.
#' @param stutter_model. Optionally a stutter_model object. See \link{global_stutter_model}.
#' @details Defines a log normal model as described by Bright et al.
#' @export
log_normal_model <- function(template, degradation = rep(0., length(template)),
                             c2, stutter_model,
                             size_regression){

  if (!is.numeric(template)){
    stop("template should be a numeric vector")
  }
  if (any(template <= 0)){
    stop("template should be positive")
  }

  if ((!is.numeric(degradation)) || (length(degradation) != length(template)) ){
    stop("degradation should be a numeric of length ", length(template))
  }

  if (any(degradation < 0)){
    stop("degradation should be non-negative")
  }

  if (!is.numeric(c2)){
    stop("c2 should be numeric")
  }
  if (length(c2) != 1){
    stop("c2 should have length 1")
  }
  if (c2 <= 0){
    stop("c2 should be positive")
  }

  model <- list()

  parameters <- list(template = template,
                     degradation = degradation,
                     c2 = c2)

  model$parameters <- parameters
  model$size_regression <- size_regression
  model$build_expected_profile <- function(...) log_normal_model_build_expected_profile(model, ...)
  model$sample_peak_heights <- function(...) log_normal_model_sample_peak_heights(model, ...)

  if (!missing(stutter_model)){
    stop("stutter is not implemented yet for log-normal model")
    model$stutter_model <- stutter_model
  }

  class(model) <- "pg_model"

  model
}

log_normal_model_build_expected_profile <- function(model, genotypes){

  if (!inherits(model, "pg_model")){
    stop("model is not of class pg_model")
  }

  parameters <- model$parameters
  size_regression <- model$size_regression

  number_of_contributors <- length(parameters$template)

  if (number_of_contributors != length(genotypes)){
    stop(number_of_contributors, " template provided for ",
         length(genotypes), " genotypes")
  }

  degradation <- parameters$degradation

  x <- data.frame(
    Marker=character(),
    Allele=character(),
    Height=numeric(),
    Size=numeric(),
    Expected=numeric(),
    stringsAsFactors = FALSE)

  for (i_contributor in seq_len(number_of_contributors)){

    template_contributor <- parameters$template[i_contributor]

    g <- genotypes[[i_contributor]]

    for (i_row in seq_len(nrow(g))){

      locus <- g$Locus[i_row]
      ab <- c(g$Allele1[i_row], g$Allele2[i_row])

      for (a in ab){
        size <- size_regression(locus, a)

        deg <- exp(-degradation[i_contributor] * (size - 100.))

        x <- add_expected_allelic_peak_height(x, locus, a, size,
                                              deg * template_contributor)
      }

    }
  }

  if (!is.null(model$stutter_model)){
    stutter_model <- model$stutter_model

    x <- stutter_model$add_expected_stutter(x)
    x$Expected <- x$ExpectedAllele + x$ExpectedStutter
  }
  else{
    x$Expected <- x$ExpectedAllele
  }

  x
}

log_normal_model_sample_peak_heights <- function(model, x){

  parameters <- model$parameters
  c2 <- parameters$c2

  x$c2 <- c2

  b <- 1000.
  x$VarianceAllele <- c2 / (b / x$ExpectedAllele + x$ExpectedAllele)

  x$HeightAllele <- 10^(log10(x$ExpectedAllele) + rnorm(n = nrow(x),mean = 0,sd = sqrt(x$VarianceAllele)))

  x$Height <- x$HeightAllele

  x
}
