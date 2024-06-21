#' @title Stutter model where the expected stutter rate depends on the allele and locus
#'
#' @param stutter_types List. See \link{stutter_type}.
#' @param size_regression Function, see \link{read_size_regression}.
#' @param sex_locus_name Character vector, defaults to "AMEL".
#' @details When a pg_model is constructed (see \link{gamma_model}), a stutter model can optionally be applied. The allele specific stutter model is commonly used with a log normal model. The expected stutter ratio for a parent allele at a locus is obtained from a linear regression of observed stutter ratios against allele length. For some loci or alleles the linear model may not be satisfactory. To override the expected stutter rates for specific alleles, a list of exceptions can be used. See \link{stutter_type} for more detail.
#' @return Object of class \code{stutter_model} to be used by e.g. \link{log_normal_model}.
#' @seealso \link{global_stutter_model} for a stutter model where the expected stutter ratio does not depend on the locus or parent allele.
#' @examples
#' # we will define an allele specific stutter model for back stutter only
#'
#' # prepare stutter regression
#' filename_bs_regression <- system.file("extdata",
#' "GlobalFiler_Stutter_3500.txt",package = "simDNAmixtures")
#' bs_regression <- read_stutter_regression(filename_bs_regression)
#'
#' # prepare exceptions, i.e. where does the regression not apply?
#' filename_bs_exceptions <- system.file("extdata",
#' "GlobalFiler_Stutter_Exceptions_3500.csv",package = "simDNAmixtures")
#' bs_exceptions <- read_stutter_exceptions(filename_bs_exceptions)
#'
#' # prepare a stutter type
#' backstutter <- stutter_type(name = "BackStutter", delta = -1,
#'                             stutter_regression = bs_regression,
#'                             stutter_exceptions = bs_exceptions)
#'
#' # assign stutter model
#' size_regression <- read_size_regression(system.file("extdata",
#' "GlobalFiler_SizeRegression.csv",package = "simDNAmixtures"))
#' bs_model <- allele_specific_stutter_model(list(backstutter), size_regression)
#' bs_model
#' @export
allele_specific_stutter_model <- function(stutter_types, size_regression, sex_locus_name = "AMEL"){

  stutter_model <- list(stutter_types = stutter_types,
                        size_regression = size_regression,
                        sex_locus_name = sex_locus_name)
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

      is_sex_locus <- marker == stutter_model$sex_locus_name
      applies_to_locus <- (!is_sex_locus) &
        (stutter$applies_to_all_loci | (marker %in% stutter$applies_to_loci))

      if (applies_to_locus){
        parent <- x$Allele[i_row]
        parent_size <- x$Size[i_row]

        if (length(stutter$delta) == 1){
          target <- get_stutter_target(parent, stutter$delta)
        }
        else{
          repeat_length <- stutter$repeat_length_by_marker[marker]
          target <- get_stutter_target(parent, stutter$delta, repeat_length)
        }

        stutter_ratio <- stutter$get_expected_stutter_ratio(marker, parent)

        # set expected stutter ratio
        x <- set_or_add_df_variable(x, marker, parent, parent_size, stutter_ratio, sr_column_name)

        # set stutter product
        x <- set_or_add_df_variable(x, marker, parent, parent_size, target, stutter_product_column_name)
      }
    }

    x$SumOfStutterRatios <- x$SumOfStutterRatios + x[[sr_column_name]]
  }

  # determine ExpectedAllele from Total Allelic Product
  x$ExpectedAllele <- x$ExpectedAllelePreStutter / (1 + x$SumOfStutterRatios)

  # determine Expected Stutters similarly
  # note that in the Bright et al. (log normal) model this is not used downstream;
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
      marker <- x_pre_stutter_products$Marker[i_row]

      if (stutter$applies_to_all_loci | (marker %in% stutter$applies_to_loci)){
        parent <- x_pre_stutter_products$Allele[i_row]

        parent_height <- x_pre_stutter_products$ExpectedAllele[i_row]

        sr <- x_pre_stutter_products[[sr_column_name]][i_row]
        target <- x_pre_stutter_products[[stutter_product_column_name]][i_row]
        target_size <- size_regression(marker, target)

        expected_stutter <- parent_height * sr / (1 + x_pre_stutter_products$SumOfStutterRatios[i_row])

        # set expected stutter
        x <- set_or_add_df_variable(x, marker, target, target_size,
                                                expected_stutter, expected_column_name)

        # store the parent
        x <- set_or_add_df_variable(x, marker, target, target_size,
                                                parent, stutter_parent_column_name)
      }
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
