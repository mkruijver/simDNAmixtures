#' @title Global stutter model where the expected stutter rate is constant across alleles and loci
#'
#' @param back_stutter_rate Numeric. (Optional)
#' @param forward_stutter_rate Numeric. (Optional)
#' @param size_regression. Function, see \link{read_size_regression}.
#' @details When a pg_model is constructed (see \link{gamma_model}), a stutter model can optionally be applied.
#' @export
global_stutter_model <- function(back_stutter_rate, forward_stutter_rate, size_regression){

  stutter_model <- list(stutters = list())
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

  x$ExpectedAllelicPreStutter <- x$ExpectedAllelic

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
      parent <- x_pre_stutter$Allele[i_row]
      target <- get_stutter_target(parent, stutter$delta)
      target_size <- size_regression(marker, target)
      parent_height <- x_pre_stutter$ExpectedAllelic[i_row]
      expected <- parent_height * stutter_rate

      # subtract stutter from Allelic
      x <- add_expected_peak_height(x, marker, parent, NA, -expected, "ExpectedAllelic")
      # add expected stutter
      x <- add_expected_peak_height(x, marker, target, target_size, expected, column_name)
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

get_stutter_target <- function(parent, delta){
  as.character(as.numeric(parent) + delta)
}

add_expected_peak_height <- function(x, marker, allele, size, expected, column_name){
  idx <- get_allele_index(x, marker, allele)

  if(length(idx)==0){

    new_df <- data.frame(Marker=marker, Allele=allele,
                         ExpectedAllelicPreStutter = 0.,
                         ExpectedAllelic = 0., ExpectedStutter = 0.,
               Size=size, stringsAsFactors = FALSE)
    new_df[[column_name]] <- expected

    return(dplyr::bind_rows(x, new_df))
  }else if(length(idx)==1){
    x[[column_name]][idx] <- x[[column_name]][idx] + expected
    return(x)
  }else{
    stop("something wrong")
  }
}

