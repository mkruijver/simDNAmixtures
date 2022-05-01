#' @title Defines a log normal model for peak height variability
#'
#' @param template Numeric vector
#' @param degradation Numeric vector of same length as template. Degradation parameters for each contributor.
#' @param LSAE Numeric vector (named) with Locus Specific Amplification Efficiencies. See \link{sample_LSAE}. Defaults to 1 for each locus.
#' @param c2 Numeric. Allele variance parameter.
#' @param k2 Optionally a numeric vector with stutter variance parameters. See \link{sample_log_normal_stutter_variance}.
#' @param model_settings List. Possible parameters: \itemize{
#'  \item locus_names. Character vector.
#'  \item degradation_parameter_cap. Numeric.
#'  \item c2_prior. Numeric of length two with shape and scale.
#'  \item LSAE_variance_prior. Numeric of length one.
#'  \item detection_threshold. Numeric vector (named) with Detection Thresholds. Defaults to 50 for each locus.
#'  \item size_regression. Function, see \link{read_size_regression}.
#'  \item stutter_model. Optionally a stutter_model object that gives expected stutter heights. See \link{global_stutter_model}.
#'  \item stutter_variability. Optionally peak height variability parameters for stutters. Required when stutter_model is supplied.
#'  }
#' @details Define a log normal model for peak height variability with the parametrisation as described by Bright et al. The model may then be used to sample DNA profiles using the \link{sample_mixture_from_genotypes} function. Alternatively, to sample many models and profiles in one go with parameters according to a specified distribution, the \link{sample_mixtures} function can be used.
#' @return Object of class \code{pg_model}.
#' @seealso \link{gamma_model}.
#' @references
#' Bright, J.A. et al. (2016). Developmental validation of STRmix™, expert software for the interpretation of forensic DNA profiles. Forensic Science International: Genetics, 23, 226-239. \doi{10.1016/j.fsigen.2016.05.007}
#' @examples
#' gf <- get_GlobalFiler_3500_data()
#' freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc.csv",
#'                            package = "simDNAmixtures"))
#'
#' k2 <- sample_log_normal_stutter_variance(gf$log_normal_settings$stutter_variability)
#'
#' model <- log_normal_model(template = 1e3, c2 = 15, k2 = k2,
#'                           model_settings = gf$log_normal_settings)
#' model
#' @export
log_normal_model <- function(template, degradation = rep(0., length(template)),
                             LSAE = stats::setNames(rep(1., length(model_settings$locus_names)),
                                             model_settings$locus_names),
                             c2, k2,
                             model_settings){

  if (!is.numeric(template)){
    stop("template should be a numeric vector")
  }
  if (any(template <= 0)){
    stop("template should be positive")
  }

  if ((!is.numeric(degradation)) || (length(degradation) != length(template)) ){
    stop("degradation should be a numeric of length ", length(template))
  }

  if (any(degradation < 0)){
    stop("degradation should be non-negative")
  }

  if (!is.numeric(c2)){
    stop("c2 should be numeric")
  }
  if (length(c2) != 1){
    stop("c2 should have length 1")
  }
  if (c2 <= 0){
    stop("c2 should be positive")
  }

  validate_log_normal_model_settings(model_settings, LSAE, k2)

  model <- list()

  parameters <- list(template = template,
                     degradation = degradation,
                     c2 = c2)

  model$locus_names <- model_settings$locus_names

  model$detection_threshold <- model_settings$detection_threshold

  model$parameters <- parameters
  model$size_regression <- model_settings$size_regression

  model$build_expected_profile_and_sample_peak_heights <- function(genotypes){
    expected_profile <- log_normal_model_build_expected_profile(model, genotypes)

    x <- log_normal_model_sample_peak_heights(model,
                                              expected_profile,
                                              model$stutter_variability)

    x$LSAE <- LSAE[x$Marker]

    x
  }

  model$sample_name_suffix <- log_normal_get_sample_name_suffix(parameters)

  if (!is.null(model_settings$stutter_model)){
    model$stutter_model <- model_settings$stutter_model
    model$stutter_variability <- model_settings$stutter_variability
    model$parameters$k2 <- k2
  }

  model$parameters$LSAE <- LSAE

  class(model) <- "pg_model"

  model
}

log_normal_model_build_expected_profile <- function(model, genotypes){

  if (!inherits(model, "pg_model")){
    stop("model is not of class pg_model")
  }

  parameters <- model$parameters
  size_regression <- model$size_regression

  number_of_contributors <- length(parameters$template)

  if (number_of_contributors != length(genotypes)){
    stop(number_of_contributors, " template provided for ",
         length(genotypes), " genotypes")
  }

  # determine deg_starts_at as minimal size
  genotypes_bound <- do.call(rbind, genotypes)
  genotypes_bound$Size1 <- genotypes_bound$Size2 <-
    numeric(nrow(genotypes_bound))

  for (i_row in seq_len(nrow(genotypes_bound))){
    genotypes_bound$Size1[i_row] <- size_regression(
      genotypes_bound$Locus[i_row], genotypes_bound$Allele1[i_row])
    genotypes_bound$Size2[i_row] <- size_regression(
      genotypes_bound$Locus[i_row], genotypes_bound$Allele2[i_row])
  }
  min_size <- min(min(genotypes_bound$Size1), min(genotypes_bound$Size2))

  degradation <- parameters$degradation

  x <- data.frame(
    Marker=character(),
    Allele=character(),
    Size=numeric(),
    Height=numeric(),
    Expected=numeric(),
    LSAE=numeric(),
    stringsAsFactors = FALSE)

  for (i_contributor in seq_len(number_of_contributors)){

    template_contributor <- parameters$template[i_contributor]

    g <- genotypes[[i_contributor]]

    for (i_row in seq_len(nrow(g))){

      locus <- g$Locus[i_row]
      ab <- c(g$Allele1[i_row], g$Allele2[i_row])

      lsae <- as.numeric(parameters$LSAE[locus])

      for (a in ab){
        size <- size_regression(locus, a)

        deg <- exp(-degradation[i_contributor] * (size - min_size))

        amount <- lsae * deg * template_contributor

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

# stutter_variability <- model$stutter_variability
# stutter_model <- model$stutter_model

log_normal_model_sample_peak_heights <- function(model, x, stutter_variability){

  parameters <- model$parameters
  size_regression <- model$size_regression

  stutter_model <- model$stutter_model

  # sample allele component of peak
  c2 <- parameters$c2

  x$c2 <- c2
  b <- 1000.
  x$VarianceAllele <- c2 / (b / x$ExpectedAllele + x$ExpectedAllele)

  x$HeightAllele <- 10^(log10(x$ExpectedAllele) + stats::rnorm(n = nrow(x),
                                                        mean = 0,
                                                        sd = sqrt(x$VarianceAllele)))

  i_stutter = 1

  # determine expected stutters based on _realised_ HeightAllele
  k2 <- parameters$k2

  for (i_stutter in seq_along(stutter_model$stutter_types)){
    stutter <- stutter_model$stutter_types[[i_stutter]]
    stutter_name <- stutter$name

    sr_max <- stutter_variability[[stutter_name]]$max_stutter_ratio

    sr_column_name <- paste0("StutterRatio", stutter_name)
    stutter_product_column_name <- paste0("StutterProduct", stutter_name)

    expected_column_name <- paste0("Expected", stutter_name)
    stutter_cap_column_name <- paste0("StutterCap", stutter_name)

    x[[stutter_cap_column_name]] <- 0.

    idx_parents <- !is.na(x[[stutter_product_column_name]])
    stutter_products <- x[[stutter_product_column_name]][idx_parents]

    idx_targets <- match(paste0(x$Marker[idx_parents], stutter_products),
                         paste0(x$Marker, x$Allele))

    x[[expected_column_name]][idx_targets] <- x$HeightAllele[idx_parents] *
                                    x[[sr_column_name]][idx_parents]

    x[[stutter_cap_column_name]][idx_targets] <-
        floor(round(x$HeightAllele[idx_parents]) * sr_max)
  }

  # sample stutter heights
  stutter_variability <- model$stutter_variability

  stutter_types <- model$stutter_model$stutter_types

  for (stutter_name in names(stutter_types)){
    expected_column <- paste0("Expected", stutter_name)
    variance_column <- paste0("Variance", stutter_name)

    height_uncapped_column <- paste0("HeightUncapped", stutter_name)

    x[[height_uncapped_column]] <- 0.

    stutter_parent_column <- paste0("StutterParent", stutter_name)

    # if inversely proportional to parent
    # then log10 O/E ~ N(0, sd = sqrt(k2 / (b / O_parent + O_parent))
    # else                            k2 / (b / E + E)

    stutter_k2 <- k2[[paste0("k2", stutter_name)]]

    if (stutter_variability[[stutter_name]]$inversely_proportional_to_parent){
      idx_stutter <- x[[expected_column]] > 0

      stutter_parent_column[idx_stutter]

      idx_parents <- match(paste0(x$Marker[idx_stutter], x[[stutter_parent_column]][idx_stutter]),
                           paste0(x$Marker, x$Allele))




      observed_parent <- x$HeightAllele[idx_parents]

      x[[variance_column]] <- 0.
      x[[variance_column]][idx_stutter] <- stutter_k2 / (b / observed_parent + observed_parent)

      x[[height_uncapped_column]][idx_stutter] <- 10^(log10(x[[expected_column]][idx_stutter]) +
                                               stats::rnorm(n = sum(idx_stutter),
                                                     mean = 0,
                                                     sd = sqrt(x[[variance_column]][idx_stutter])))

    }
    else{
      x[[variance_column]] <- stutter_k2 / (b / x[[expected_column]] + x[[expected_column]])

      x[[height_uncapped_column]] <- 10^(log10(x[[expected_column]]) + stats::rnorm(n = nrow(x),
                                                     mean = 0,
                                                     sd = sqrt(x[[variance_column]])))
    }


  }

  if (!is.null(stutter_types)){
    # add up stutters
    x$HeightUncappedStutter <- rowSums(x[paste0("HeightUncapped", names(stutter_types))])

    x$StutterCap <- do.call(pmax, x[paste0("StutterCap", names(stutter_types))])

    x$HeightStutter <- pmin(x$HeightUncappedStutter, x$StutterCap)

    x$Height <- x$HeightAllele + x$HeightStutter

    # fix expected stutter total
    x$ExpectedStutter <- rowSums(x[paste0("Expected", names(stutter_types))])
  }
  else{
    x$Height <- x$HeightAllele
    x$ExpectedStutter <- 0.
  }

  x$Expected <- x$ExpectedStutter + x$ExpectedAllele

  # add detection threshold
  x$DetectionThreshold <- model$detection_threshold[x$Marker]
  x$HeightAtOrAboveDetectionThreshold <- round(x$Height) >= x$DetectionThreshold

  x
}

log_normal_get_sample_name_suffix <- function(parameters){

  noc_label <- length(parameters$template)
  template_label <- paste(round(parameters$template,digits = 0), collapse = "_")

  number_of_buckets <- 5
  buckets <- findInterval(x = parameters$degradation/0.01,
                          vec = seq(from = 0, to = 1,
                                    length = number_of_buckets + 1))
  deg_label <- paste0(letters[buckets],collapse = "")

  paste0("N", noc_label, "_", template_label,"_",deg_label)
}
