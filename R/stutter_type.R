#' @title Defines a stutter type to be used in the allele specific stutter model.
#'
#' @param name Character. Name of the stutter, e.g. "BackStutter"
#' @param delta Numeric. When length one, repeat units gained (lost when negative). When length two, the second element is the number of base pairs gained (lost).
#' @param applies_to_all_loci Logical. Defaults to TRUE.
#' @param stutter_regression Function. See \link{read_stutter_regression}.
#' @param stutter_exceptions Optionally a list. See \link{read_stutter_exceptions}.
#' @param applies_to_loci Optionally a character vector of locus names to which this stutter type applies.
#' @param repeat_length_by_marker Optionally a named integer vector with repeat lengths by marker. Only needed when delta is of length two.
#' @details When a pg_model is constructed (see \link{log_normal_model}), a stutter model can optionally be applied.
#' @return Object of class \code{stutter_type} to be passed to \link{allele_specific_stutter_model}.
#' @examples
#' filename_bs_exceptions <- system.file("extdata",
#' "GlobalFiler_Stutter_Exceptions_3500.csv",package = "simDNAmixtures")
#' bs_exceptions <- read_stutter_exceptions(filename_bs_exceptions)
#'
#' filename_bs_regression <- system.file("extdata",
#' "GlobalFiler_Stutter_3500.txt",package = "simDNAmixtures")
#' bs_regression <- read_stutter_regression(filename_bs_regression)
#'
#' backstutter <- stutter_type(name = "BackStutter", delta = -1,
#'                             stutter_regression = bs_regression,
#'                             stutter_exceptions = bs_exceptions)
#' @export
stutter_type <- function(name, delta,
                         applies_to_all_loci = TRUE,
                         stutter_regression,
                         stutter_exceptions,
                         applies_to_loci,
                         repeat_length_by_marker){

  stutter <- list()
  class(stutter) <- "stutter_type"

  if (!is.character(name)){
    stop("name is not a character")
  }
  if (length(name)!=1){
    stop("name is not length 1")
  }

  if ((!is.logical(applies_to_all_loci)) || (length(applies_to_all_loci) != 1)){
    stop("applies_to_all_loci needs to be a logical of length 1")
  }

  if ((!applies_to_all_loci) & missing(applies_to_loci)){
    stop("applies_to_all_loci is FALSE but applies_to_loci is missing")
  }

  if (!missing(applies_to_loci)){
    if (!is.character(applies_to_loci)){
      stop("applies_to_loci is not a character vector of locus names")
    }
  }

  if (!is.numeric(delta)){
    stop("delta is not numeric")
  }

  if ((length(delta) != 1) & (length(delta) != 2)){
    stop("delta is not length 1 or 2")
  }

  if ((length(delta) == 2) & missing(repeat_length_by_marker)){
    stop("repeat_length_by_marker is missing and is needed because delta is length 2")
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
  stutter$applies_to_all_loci <- applies_to_all_loci

  if (!missing(applies_to_loci)){
    stutter$applies_to_loci <- applies_to_loci
  }

  if (!missing(repeat_length_by_marker)){
    stutter$repeat_length_by_marker <- repeat_length_by_marker
  }

  stutter$regression <- stutter_regression

  if (!missing(stutter_exceptions)){
    stutter$exceptions <- stutter_exceptions
  }

  stutter$get_expected_stutter_ratio <- function(locus, allele){

    exception <- NULL
    if (!is.null(stutter$exceptions[[locus]])){
      if (as.character(allele) %in% names(stutter$exceptions[[locus]])){
        exception <- stutter$exceptions[[locus]][[as.character(allele)]]
      }
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
