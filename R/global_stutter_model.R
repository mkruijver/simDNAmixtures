#' @title Global stutter model where the expected stutter rate is constant across alleles and loci
#'
#' @param back_stutter_rate Numeric. (Optional)
#' @param forward_stutter_rate Numeric. (Optional)
#' @param size_regression Function, see \link{read_size_regression}.
#' @param sex_locus_name Character vector, defaults to "AMEL".
#' @details When a pg_model is constructed (see \link{gamma_model}), a stutter model can optionally be applied. In the global stutter model, the expected stutter rate is constant across all loci and for all parent alleles.
#' @return Object of class \code{stutter_model} to be used by e.g. \link{gamma_model}.
#' @seealso \link{allele_specific_stutter_model} for a stutter model where the expected stutter rate depends on the allele and locus.
#' @examples
#' # the stutter model needs a size regression to determine fragment length
#' # of stutter products
#' size_regression <- read_size_regression(system.file("extdata",
#' "GlobalFiler_SizeRegression.csv",package = "simDNAmixtures"))
#'
#' # define a stutter model with an expected back stutter rate of 10%
#' stutter_model <- global_stutter_model(back_stutter_rate = 0.1,
#'                                      size_regression = size_regression)
#'
#' stutter_model
#' @export
global_stutter_model <- function(back_stutter_rate, forward_stutter_rate,
                                 size_regression, sex_locus_name = "AMEL"){

  stutter_model <- list(stutters = list(),
                        sex_locus_name = sex_locus_name)
  class(stutter_model) <- "stutter_model"

  if (!missing(back_stutter_rate)){
    if (!(is.numeric(back_stutter_rate) & (length(back_stutter_rate)==1) )){
      stop("back_stutter_rate needs to be a numeric of length 1")
    }
    if (back_stutter_rate < 0 | back_stutter_rate > 1){
      stop("back_stutter_rate needs to be between 0 and 1")
    }

    stutter_model$stutters$BackStutter <- list(
      rate = back_stutter_rate,
      delta = -1.
    )
  }

  if (!missing(forward_stutter_rate)){
    if (!(is.numeric(forward_stutter_rate) & (length(forward_stutter_rate)==1) )){
      stop("forward_stutter_rate needs to be a numeric of length 1")
    }
    if (forward_stutter_rate < 0 | forward_stutter_rate > 1){
      stop("forward_stutter_rate needs to be between 0 and 1")
    }

    stutter_model$stutters$ForwardStutter <- list(
      rate = forward_stutter_rate,
      delta = 1.
    )
  }

  stutter_model$size_regression <- size_regression
  stutter_model$add_expected_stutter <- function(...) global_stutter_model_add_expected_stutter(stutter_model, ...)

  stutter_model
}

global_stutter_model_add_expected_stutter <- function(stutter_model, x){

  x_pre_stutter <- x

  x$ExpectedAllelePreStutter <- x$ExpectedAllele

  x$ExpectedStutter <- 0.

  size_regression <- stutter_model$size_regression

  # apply all stutter types in the model
  for (i_stutter in seq_along(stutter_model$stutters)){
    stutter <- stutter_model$stutters[[i_stutter]]

    stutter_name <- names(stutter_model$stutters)[[i_stutter]]

    stutter_rate <- stutter$rate
    column_name <- paste0("Expected", stutter_name)
    x[[column_name]] <- 0.

    for (i_row in seq_len(nrow(x_pre_stutter))){
      marker <- x_pre_stutter$Marker[i_row]

      applies_to_locus <- marker != stutter_model$sex_locus_name

      if (applies_to_locus){
        parent <- x_pre_stutter$Allele[i_row]
        target <- get_stutter_target(parent, stutter$delta)
        target_size <- size_regression(marker, target)
        parent_height <- x_pre_stutter$ExpectedAllele[i_row]
        expected <- parent_height * stutter_rate

        # subtract stutter from Allele
        x <- add_expected_peak_height(x, marker, parent, NA, -expected, "ExpectedAllele")
        # add expected stutter
        x <- add_expected_peak_height(x, marker, target, target_size, expected, column_name)
      }
    }

    x$ExpectedStutter <- x$ExpectedStutter + x[[column_name]]
  }

  # set NAs from rowbinding to zero
  for (stutter_name in names(stutter_model$stutters)){
    column_name <- paste0("Expected", stutter_name)

    x[[column_name]][is.na(x[[column_name]])] <- 0.
  }

  x
}
