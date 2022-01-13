#' @title Defines a stutter type to be used in an allele specific stutter model.
#'
#' @param name Character. Name of the stutter, e.g. "BackStutter"
#' @param delta Numeric. Repeat units gained (lost when negative).
#' @param stutter_regression Function. See \link{read_stutter_regression}.
#' @param stutter_exceptions Optionally a list. See \link{read_stutter_exceptions}.
#' @details When a pg_model is constructed (see \link{gamma_model}), a stutter model can optionally be applied.
#' @examples
#' filename_bs_exceptions <- system.file("extdata","GlobalFiler_Stutter_Exceptions_3500.csv",package = "SimMixDNA")
#' bs_exceptions <- read_stutter_exceptions(filename_bs_exceptions)
#'
#' filename_bs_regression <- system.file("extdata","GlobalFiler_Stutter_3500.txt",package = "SimMixDNA")
#' bs_regression <- read_stutter_regression(filename_bs_regression)
#' backstutter <- stutter_type(name = "BackStutter", delta = -1,
#'                             stutter_regression = bs_regression,
#'                             stutter_exceptions = bs_exceptions)
#' @export
stutter_type <- function(name, delta, stutter_regression, stutter_exceptions){

  stutter <- list()
  class(stutter) <- "stutter_type"

  if (!is.character(name)){
    stop("name is not a character")
  }
  if (length(name)!=1){
    stop("name is not length 1")
  }
  if (!is.numeric(delta)){
    stop("delta is not numeric")
  }
  if (length(delta) != 1){
    stop("delta is not length 1")
  }
  if (!is.function(stutter_regression)){
    stop("stutter_regression is not a function")
  }
  if (!missing(stutter_exceptions)){
    if (!is.list(stutter_exceptions)){
      stop("stutter_exceptions is not a list")
    }
  }

  stutter$name <- name
  stutter$delta <- delta

  stutter$regression <- stutter_regression

  if (!missing(stutter_exceptions)){
    stutter$exceptions <- stutter_exceptions
  }

  stutter$get_expected_stutter_ratio <- function(locus, allele){

    if (!is.null(stutter$exceptions)){
      exception <- stutter$exceptions[[locus]][[as.character(allele)]]
    }

    if (isTRUE(exception > 0)){
      return(exception)
    }
    else{
      return(stutter$regression(locus, allele))
    }
  }

  stutter
}

