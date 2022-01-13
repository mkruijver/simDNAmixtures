#' @title Stutter model where the expected stutter rate depends on the allele and locus
#'
#' @param stutter_types List. See \link{stutter_type}.
#' @param size_regression. Function, see \link{read_size_regression}.
#' @details When a pg_model is constructed (see \link{gamma_model}), a stutter model can optionally be applied.
#' @export
allele_specific_stutter_model <- function(stutter_types, size_regression){

  stutter_model <- list(stutter_types = stutter_types,
                        size_regression = size_regression)
  class(stutter_model) <- "stutter_model"

  stutter_model$add_expected_stutter <- function(...) allele_specific_stutter_model_add_expected_stutter(stutter_model, ...)

  stutter_model
}

allele_specific_stutter_model_add_expected_stutter <- function(stutter_model, x){

  x_pre_stutter <- x

  x$ExpectedAllelePreStutter <- x$ExpectedAllele

  x$ExpectedStutter <- 0.

  size_regression <- stutter_model$size_regression

  # compute ExpectedStutter for all stutter types in the model
  for (i_stutter in seq_along(stutter_model$stutter_types)){
    stutter <- stutter_model$stutter_types[[i_stutter]]

    stutter_name <- stutter$name

    column_name <- paste0("Expected", stutter_name)
    x[[column_name]] <- 0.

    for (i_row in seq_len(nrow(x_pre_stutter))){
      marker <- x_pre_stutter$Marker[i_row]
      parent <- x_pre_stutter$Allele[i_row]
      target <- get_stutter_target(parent, stutter$delta)
      target_size <- size_regression(marker, target)

      parent_height <- x_pre_stutter$ExpectedAllele[i_row]

      stutter_rate <- stutter$get_expected_stutter_ratio(marker, parent)
      expected <- parent_height * stutter_rate

      # subtract stutter from Allele
      x <- add_expected_peak_height(x, marker, parent, NA, -expected, "ExpectedAllele")
      # add expected stutter
      x <- add_expected_peak_height(x, marker, target, target_size, expected, column_name)
    }

    x$ExpectedStutter <- x$ExpectedStutter + x[[column_name]]
  }

  # set NAs from rowbinding to zero
  for (stutter_name in names(stutter_model$stutter_types)){
    column_name <- paste0("Expected", stutter_name)

    x[[column_name]][is.na(x[[column_name]])] <- 0.
  }

  x
}
