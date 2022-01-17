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

  x$ExpectedAllelePreStutter <- x$ExpectedAllele

  x$Expected <- x$ExpectedAllele <- x$ExpectedStutter <- x$SumOfStutterRatios <- 0.

  size_regression <- stutter_model$size_regression

  # compute expected stutter rates and determine stutter products
  for (i_stutter in seq_along(stutter_model$stutter_types)){
    stutter <- stutter_model$stutter_types[[i_stutter]]

    stutter_name <- stutter$name

    sr_column_name <- paste0("StutterRatio", stutter_name)
    x[[sr_column_name]] <- 0.

    stutter_product_column_name <- paste0("StutterProduct", stutter_name)

    stutter_parent_column_name <- paste0("StutterParent", stutter_name)
    x[[stutter_product_column_name]] <- x[[stutter_parent_column_name]] <- NA_character_

    for (i_row in seq_len(nrow(x))){
      marker <- x$Marker[i_row]
      parent <- x$Allele[i_row]
      parent_size <- x$Size[i_row]

      target <- SimMixDNA:::get_stutter_target(parent, stutter$delta)
      stutter_ratio <- stutter$get_expected_stutter_ratio(marker, parent)

      # set expected stutter ratio
      x <- SimMixDNA:::set_or_add_df_variable(x, marker, parent, parent_size, stutter_ratio, sr_column_name)

      # set stutter product
      x <- SimMixDNA:::set_or_add_df_variable(x, marker, parent, parent_size, target, stutter_product_column_name)
    }

    x$SumOfStutterRatios <- x$SumOfStutterRatios + x[[sr_column_name]]
  }

  # determine ExpectedAllele from Total Allelic Product
  x$ExpectedAllele <- x$ExpectedAllelePreStutter / (1 + x$SumOfStutterRatios)

  # determine Expected Stutters similarly
  # note that in the Bright et al. model this is not used downstream;
  # rather the Expected Stutter is overwritten as being proportional
  # to the _observed_ parent height which is not available at this stge
  x_pre_stutter_products <- x # backup and grow x in place

  expected_stutter_column_names <- character()

  i_stutter=1
  for (i_stutter in seq_along(stutter_model$stutter_types)){
    stutter <- stutter_model$stutter_types[[i_stutter]]

    stutter_name <- stutter$name

    sr_column_name <- paste0("StutterRatio", stutter_name)
    stutter_product_column_name <- paste0("StutterProduct", stutter_name)
    stutter_parent_column_name <- paste0("StutterParent", stutter_name)

    expected_column_name <- paste0("Expected", stutter_name)
    x[[expected_column_name]] <- 0.
    expected_stutter_column_names <- c(expected_stutter_column_names, expected_column_name)

    i_row=1
    for (i_row in seq_len(nrow(x_pre_stutter_products))){
      marker <- x$Marker[i_row]
      parent <- x$Allele[i_row]

      parent_height <- x$ExpectedAllele[i_row]

      sr <- x[[sr_column_name]][i_row]
      target <- x[[stutter_product_column_name]][i_row]
      target_size <- size_regression(marker, target)

      expected_stutter <- parent_height * sr / (1 + x$SumOfStutterRatios[i_row])

      # set expected stutter
      x <- SimMixDNA:::set_or_add_df_variable(x, marker, target, target_size,
                                              expected_stutter, expected_column_name)

      # store the parent
      x <- SimMixDNA:::set_or_add_df_variable(x, marker, target, target_size,
                                              parent, stutter_parent_column_name)

    }
  }

  # set NAs from rowbinding to zero
  for (column_name in names(x)[sapply(x, is.numeric)]){
    x[[column_name]][is.na(x[[column_name]])] <- 0.
  }

  x$ExpectedStutter <- rowSums(x[expected_stutter_column_names])

  x$Expected <- x$ExpectedStutter + x$ExpectedAllele

  x
}
