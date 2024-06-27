#' @title Defines a (semi-continuous) drop model
#'
#' @param dropout_probabilities Numeric vector with values between 0 and 1. Dropout probabilities for each contributor.
#' @param drop_in_rate Numeric vector of length one. Expected number of drop-ins per locus. Default is 0.
#' @param freqs Optionally a list with allele frequencies (needed when drop_in_rate > 0). See \link{read_allele_freqs}.
#' @param model_settings List. Possible parameters: \itemize{
#'  \item locus_names. Character vector.
#'  \item size_regression. Function, see \link{read_size_regression}.
#'  }
#' @details Define the classes semi-continuous drop-model. The model may then be used to sample DNA profiles using the \link{sample_mixture_from_genotypes} function. Alternatively, to sample many models and profiles in one go with parameters according to a specified distribution, the \link{sample_mixtures} function can be used.
#' @return Object of class \code{pg_model}.
#' @seealso \link{gamma_model}, \link{log_normal_model}.
#' @references
#' Slooten, K. (2017). Accurate assessment of the weight of evidence for DNA mixtures by integrating the likelihood ratio. Forensic Science International: Genetics, 27, 1-16. \doi{10.1016/j.fsigen.2016.11.001}
#' @examples
#' gf <- gf_configuration()
#' freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc_022024.csv",
#'                            package = "simDNAmixtures"))
#'
#' k2 <- sample_log_normal_stutter_variance(gf$log_normal_settings$stutter_variability)
#'
#' model <- log_normal_model(template = 1e3, c2 = 15, k2 = k2,
#'                           model_settings = gf$log_normal_settings)
#' model
#' @export
drop_model <- function(dropout_probabilities, drop_in_rate = 0.,
                             freqs,
                             model_settings){

  .validate_clamp(dropout_probabilities,
                  minimum = 0., maximum = 1.,
                  strict_maximum = TRUE)
  .validate_numeric(drop_in_rate, require_nonnegative = TRUE)

  model <- list()

  parameters <- list(model = "drop_model",
                     dropout_probabilities = dropout_probabilities,
                     drop_in_rate = drop_in_rate)

  if (!missing(freqs)){
    parameters$freqs = freqs
  }

  model$locus_names <- model_settings$locus_names
  model$parameters <- parameters
  model$size_regression <- model_settings$size_regression

  model$build_expected_profile_and_sample_peak_heights <- function(genotypes){
    expected_profile <- drop_model_build_expected_profile(model, genotypes)

    x <- drop_model_sample_dropouts(model, expected_profile)
    x
  }

  model$sample_name_suffix <- drop_model_get_sample_name_suffix(parameters)
  class(model) <- "pg_model"

  model
}

drop_model_build_expected_profile <- function(model, genotypes){

  if (!inherits(model, "pg_model")){
    stop("model is not of class pg_model")
  }

  parameters <- model$parameters
  size_regression <- model$size_regression

  number_of_contributors <- length(parameters$dropout_probabilities)

  if (number_of_contributors != length(genotypes)){
    stop(number_of_contributors, " dropout probabilities provided for ",
         length(genotypes), " genotypes")
  }

  x <- data.frame(
    Marker=character(),
    Allele=character(),
    Size=numeric(),
    Height=numeric(),
    stringsAsFactors = FALSE)

  for (i_contributor in seq_len(number_of_contributors)){
    x[[paste0("n", i_contributor)]] <- rep(0L, nrow(x))

    # extract allele columns
    g <- genotypes[[i_contributor]]
    allele_columns <- .get_allele_columns(g)

    for (i_row in seq_len(nrow(allele_columns))){

      locus <- g$Locus[i_row]

      for (i_allele in seq_len(ncol(allele_columns))){
        a <- allele_columns[i_row, i_allele]

        if (!is.na(a)){
          size <- size_regression(locus, a)

          x <- set_or_add_df_variable(x, locus, a, size,
                                      1L, paste0("n", i_contributor),
                                      sum = TRUE)
        }
      }
    }
  }

  x
}

drop_model_sample_dropouts <- function(model, x){

  parameters <- model$parameters

  dropout_probabilities <- parameters$dropout_probabilities
  number_of_contributors <- length(dropout_probabilities)

  # sample allele component of peak
  for (i_contributor in seq_len(number_of_contributors)){
    n_col <- paste0("n", i_contributor)

    not_dropout_probability <- 1.0 - dropout_probabilities[i_contributor]

    # sample number of remaining alleles
    x[[n_col]][is.na(x[[n_col]])] <- 0L

    x[[paste0("x", i_contributor)]] <- stats::rbinom(n = nrow(x), size = x[[n_col]],
                                              prob = not_dropout_probability)
  }

  x$n <- rowSums(x[paste0("n", seq_len(number_of_contributors))])
  x$x <- rowSums(x[paste0("x", seq_len(number_of_contributors))])

  x$Height <- ifelse(x$x > 0, yes = 1000., no = 0.)
  x$HeightAtOrAboveDetectionThreshold <- x$Height > 0.

  x
}

drop_model_get_sample_name_suffix <- function(parameters){

  noc_label <- length(parameters$dropout_probabilities)

  number_of_buckets <- 26
  buckets <- findInterval(x = parameters$dropout_probabilities,
                          vec = seq(from = 0, to = 1,
                                    length = number_of_buckets + 1))
  drop_label <- paste0(LETTERS[buckets],collapse = "")

  paste0("N", noc_label, "_", drop_label)
}
