#' @title Defines a gamma model for peak height variability
#'
#' @param mixture_proportions Numeric vector with the mixture proportion for each contributor.
#' @param mu Numeric. Expectation of a full heterozygote contributing allele peak height.
#' @param cv Numeric. Coefficient of variation of a full heterozygote contributing allele peak height
#' @param degradation_beta Numeric Vector of same length as mixture_proportions. Degradation slope parameters for each contributor. Defaults to 1 for each contributor (i.e. not degraded)
#' @param LSAE Numeric vector (named) with Locus Specific Amplification Efficiencies. See \link{sample_LSAE}. Defaults to 1 for each locus.
#' @param model_settings List. Possible parameters: \itemize{
#'  \item locus_names. Character vector.
#'  \item detection_threshold. Numeric vector (named) with Detection Thresholds.
#'  \item size_regression. Function, see \link{read_size_regression}.
#'  \item stutter_model. Optionally a stutter_model object that gives expected stutter heights. See \link{global_stutter_model}.
#'  }
#' @details Define a gamma model for peak height variability with the parametrisation as described by Bleka et al. The model may then be used to sample DNA profiles using the \link{sample_mixture_from_genotypes} function. Alternatively, to sample many models and profiles in one go with parameters according to a specified distribution, the \link{sample_mixtures} function can be used.
#' @return Object of class \code{pg_model}.
#' @seealso \link{log_normal_model}.
#' @examples
#' # read allele frequencies
#' freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc_022024.csv",
#'                                        package = "simDNAmixtures"))
#'
#' gf <- gf_configuration()
#'
#' # define the gamma model for peak heights
#' model <- gamma_model(mixture_proportions = 1, mu = 1000.,
#'                     cv = 0.1, model_settings = gf$gamma_settings_no_stutter)
#'
#' # sample a single source profile (1-person 'mixture')
#' u1 <- sample_contributor_genotypes("U1", freqs, loci = gf$autosomal_markers)
#' sample <- sample_mixture_from_genotypes(u1, model)
#'
#' # peaks follow a gamma distribution with an expected height of
#' # 1,000 for heterozygous alleles; 2,000 for homozygotes
#' hist(sample$Height)
#'
#' # the gamma distribution is more obvious if many samples are taken
#' many_samples <- replicate(n = 1e2,
#'                           sample_mixture_from_genotypes(u1, model),
#'                           simplify = FALSE)
#'
#' hist(sapply(many_samples, function(x) x$Height))
#' @references
#' Bleka, Ø., Storvik, G., & Gill, P. (2016). EuroForMix: An open source software based on a continuous model to evaluate STR DNA profiles from a mixture of contributors with artefacts. Forensic Science International: Genetics, 21, 35-44. \doi{10.1016/j.fsigen.2015.11.008}
#' @export
gamma_model <- function(mixture_proportions, mu, cv,
                        degradation_beta = rep(1., length(mixture_proportions)),
                        LSAE = stats::setNames(rep(1., length(model_settings$locus_names)),
                                        model_settings$locus_names),
                        model_settings){

  .validate_numeric(mixture_proportions,require_strictly_positive = TRUE,
                    require_length_one = FALSE)

  if (!isTRUE(all.equal(sum(mixture_proportions), 1))){
    stop("mixture_proportions should sum to 1")
  }

  if (length(degradation_beta) != length(mixture_proportions)){
    stop("degradation_beta should be a numeric of length ",
         length(mixture_proportions))
  }
  .validate_numeric(degradation_beta, require_nonnegative = TRUE,
                    require_length_one = FALSE)
  if (any(degradation_beta > 1)){
    stop("degradation_beta should not exceed 1")
  }

  locus_names <- model_settings$locus_names
  detection_threshold <- model_settings$detection_threshold
  size_regression <- model_settings$size_regression
  stutter_model <- model_settings$stutter_model

  if (!is.character(locus_names)){
    stop("locus_names needs to be a character vector")
  }

  if (!all(locus_names %in% names(LSAE))){
    stop("all locus names need to be in names(LSAE)")
  }

  if (!all(locus_names %in% names(detection_threshold))){
    stop("all locus names need to be in names(detection_threshold)")
  }

  .validate_numeric(LSAE, require_nonnegative = TRUE, require_length_one = FALSE)
  .validate_numeric(detection_threshold, require_nonnegative = TRUE,
                    require_length_one = FALSE)
  .validate_numeric(mu, require_strictly_positive = TRUE,
                    require_length_one = TRUE)
  .validate_numeric(cv, require_nonnegative = TRUE,
                    require_length_one = TRUE)

  if (!is.null(stutter_model)){
    if (!inherits(stutter_model, "stutter_model")){
      stop("stutter_model is not of class stutter_model")
    }
  }

  model <- list()

  model$locus_names <- locus_names
  model$detection_threshold <- detection_threshold

  parameters <- list(model = "gamma_model",
                     mixture_proportions = mixture_proportions,
                     mu = mu,
                     cv = cv,
                     degradation_beta = degradation_beta,
                     LSAE = LSAE)

  model$parameters <- parameters

  model$sample_name_suffix <- gamma_get_sample_name_suffix(parameters)

  model$size_regression <- size_regression

  model$build_expected_profile_and_sample_peak_heights <- function(genotypes){
    expected_profile <- gamma_model_build_expected_profile(model, genotypes)
    x <- gamma_model_sample_peak_heights(model, expected_profile)

    x$LSAE <- LSAE[x$Marker]

    x
  }

  if (!missing(stutter_model)){
    model$stutter_model <- stutter_model
  }

  class(model) <- "pg_model"

  model
}

gamma_model_build_expected_profile <- function(model, genotypes){

  if (!inherits(model, "pg_model")){
    stop("model is not of class pg_model")
  }

  parameters <- model$parameters
  size_regression <- model$size_regression

  number_of_contributors <- length(parameters$mixture_proportions)

  if (number_of_contributors != length(genotypes)){
    stop(number_of_contributors, " mixture proportions provided for ",
         length(genotypes), " genotypes")
  }

  beta <- parameters$degradation_beta

  x <- data.frame(
    Marker=character(),
    Allele=character(),
    Size=numeric(),
    Height=numeric(),
    Expected=numeric(),
    LSAE=numeric(),
    stringsAsFactors = FALSE)

  for (i_contributor in seq_len(number_of_contributors)){

    mu_contributor <- parameters$mu * parameters$mixture_proportions[i_contributor]

    g <- genotypes[[i_contributor]]

    # extract allele columns
    allele_columns <- .get_allele_columns(g)

    for (i_row in seq_len(nrow(allele_columns))){

      locus <- g$Locus[i_row]

      lsae <- as.numeric(parameters$LSAE[locus])

      for (i_allele in seq_len(ncol(allele_columns))){
        a <- allele_columns[i_row, i_allele]

        if (!is.na(a)){
          size <- size_regression(locus, a)

          deg <- beta[i_contributor] ^ ((size - 125.) / 100.)

          amount <- lsae * deg * mu_contributor
          x <- add_expected_allelic_peak_height(x, locus, a, size, amount)
        }
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

gamma_model_sample_peak_heights <- function(model, x){

  parameters <- model$parameters

  # scale is constant (mu * cv^2)
  # shape is determined by expected peak height
  # by using that for a gamma distribution,
  # the mean is equal to shape*scale
  # so shape is obtained as expected / scale

  cv <- parameters$cv

  x$cv <- cv
  x$Scale <- parameters$mu * cv^2
  x$Shape <- x$Expected / x$Scale

  x$Height <- stats::rgamma(n = nrow(x), shape = x$Shape, scale = x$Scale)

  # add detection threshold
  x$DetectionThreshold <- model$detection_threshold[x$Marker]
  x$HeightAtOrAboveDetectionThreshold <- round(x$Height) >= x$DetectionThreshold

  x
}

gamma_get_sample_name_suffix <- function(parameters){

  noc_label <- length(parameters$mixture_proportions)

  mu_label <- round(parameters$mu)

  mixture_proportions_label <- paste(round(100 * parameters$mixture_proportions),
                          collapse = "_")

  number_of_buckets <- 5
  buckets <- findInterval(x = 1-parameters$degradation,
                          vec = seq(from = 0, to = 1,
                                    length = number_of_buckets + 1))
  deg_label <- paste0(letters[buckets],collapse = "")

  paste0("N", noc_label, "_", mu_label, "_", mixture_proportions_label,"_",deg_label)
}
