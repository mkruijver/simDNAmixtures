#' @title Defines a gamma model for peak height variability
#'
#' @param mixture_proportions Numeric vector
#' @param mu Numeric. Expectation of a full heterozygote contributing allele peak height.
#' @param cv Numeric. Coefficient of variation of a full heterozygote contributing allele peak height
#' @param degradation_beta Numeric Vector of same length as mixture_proportions. Degradation slope parameters for each contributor.
#' @param LSAE Numeric vector (named) with Locus Specific Amplification Efficiencies. See \link{sample_LSAE}. Defaults to 1 for each locus.
#' @param model_settings List. Possible parameters: \itemize{
#'  \item locus_names. Character vector.
#'  \item detection_threshold. Numeric vector (named) with Detection Thresholds. Defaults to 50 for each locus.
#'  \item size_regression. Function, see \link{read_size_regression}.
#'  \item stutter_model. Optionally a stutter_model object that gives expected stutter heights. See \link{global_stutter_model}.
#'  }
#' @details Defines a gamma model as described by Bleka et al.
#' @export
gamma_model <- function(mixture_proportions, mu, cv,
                        degradation_beta = rep(1., length(mixture_proportions)),
                        LSAE = setNames(rep(1., length(model_settings$locus_names)),
                                        model_settings$locus_names),
                        model_settings){

  if (!is.numeric(mixture_proportions)){
    stop("mixture_proportions should be a numeric vector")
  }
  if (any(mixture_proportions <= 0)){
    stop("mixture_proportions should be positive")
  }
  if (!isTRUE(all.equal(sum(mixture_proportions), 1))){
    stop("mixture_proportions should sum to 1")
  }

  if ((!is.numeric(degradation_beta)) || (length(degradation_beta) != length(mixture_proportions)) ){
    stop("degradation_beta should be a numeric of length ", length(mixture_proportions))
  }

  if (any(degradation_beta <= 0)){
    stop("degradation_beta should be positive")
  }

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

  if (!is.numeric(LSAE)){
    stop("LSAE needs to be numeric")
  }

  if (!is.numeric(detection_threshold)){
    stop("detection_threhold needs to be numeric")
  }

  if ((!is.numeric(mu)) || (length(mu) != 1) ){
    stop("mu should be a numeric of length 1")
  }
  if (mu <= 0){
    stop("mu should be positive")
  }

  if ((!is.numeric(cv)) || (length(cv) != 1) ){
    stop("cv should be a numeric of length 1")
  }
  if (cv < 0){
    stop("cv should be non-negative")
  }

  if (!is.null(stutter_model)){
    if (!inherits(stutter_model, "stutter_model")){
      stop("stutter_model is not of class stutter_model")
    }
  }

  model <- list()

  model$locus_names <- locus_names
  model$LSAE <- LSAE
  model$detection_threshold <- detection_threshold

  parameters <- list(mixture_proportions = mixture_proportions,
                     mu = mu,
                     cv = cv,
                     degradation_beta = degradation_beta)

  model$parameters <- parameters
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

    for (i_row in seq_len(nrow(g))){

      locus <- g$Locus[i_row]
      ab <- c(g$Allele1[i_row], g$Allele2[i_row])

      lsae <- as.numeric(model$LSAE[locus])

      for (a in ab){
        size <- size_regression(locus, a)

        deg <- beta[i_contributor] ^ ((size - 125.) / 100.)

        amount <- lsae * deg * mu_contributor
        x <- add_expected_allelic_peak_height(x, locus, a, size, amount)
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

  # shape is constant (1 / cv^2)
  # scale is determined by expected peak height

  cv <- parameters$cv

  x$cv <- cv
  x$Shape <- 1 / (cv*cv)
  x$Scale <- cv*cv * x$Expected

  x$Height <- rgamma(n = nrow(x), shape = x$Shape, scale = x$Scale)

  # add detection threshold
  x$DetectionThreshold <- model$detection_threshold[x$Marker]
  x$HeightAtOrAboveDetectionThreshold <- round(x$Height) >= x$DetectionThreshold

  x
}
