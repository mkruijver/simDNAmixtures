#' @title Defines a log normal model for peak height variability
#'
#' @param template Numeric vector
#' @param degradation Numeric vector of same length as template. Degradation parameters for each contributor.
#' @param locus_names Character vector.
#' @param LSAE Numeric vector (named) with Locus Specific Amplification Efficiencies. See \link{sample_LSAE}. Defaults to 1 for each locus.
#' @param c2 Numeric. Allele variance parameter.
#' @param k2 Optionally a numeric vector with stutter variance parameters. See \link{sample_log_normal_stutter_variance}.
#' @param size_regression Function, see \link{read_size_regression}.
#' @param stutter_model Optionally a stutter_model object that gives expected stutter heights. See \link{global_stutter_model}.
#' @param stutter_variability Optionally peak height variability parameters for stutters. Required when stutter_model is supplied.
#' @details Defines a log normal model as described by Bright et al.
#' @export
log_normal_model <- function(template, degradation = rep(0., length(template)),
                             locus_names,
                             LSAE = setNames(rep(1., length(locus_names)), locus_names),
                             c2, k2, stutter_model,
                             stutter_variability,
                             size_regression){

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

  if (!is.character(locus_names)){
    stop("locus_names needs to be a character vector")
  }

  if (!all(locus_names %in% names(LSAE))){
    stop("all locus names need to be in names(LSAE)")
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

  if (!missing(stutter_model)){
    if (missing(stutter_variability)){
      stop("stutter_variability needs to be supplied when stutter_model is supplied")
    }

    if (missing(k2)){
      stop("k2 needs to be supplied when stutter_model is supplied")
    }

    expected_k2_names <- paste0("k2", names(stutter_model$stutter_types))

    if (!identical(expected_k2_names, names(k2))){
      stop("k2 does not have expected names: ", paste(expected_k2_names, collapse = ", "))
    }
  }

  model <- list()

  parameters <- list(template = template,
                     degradation = degradation,
                     c2 = c2)

  model$locus_names <- locus_names
  model$LSAE <- LSAE

  model$parameters <- parameters
  model$size_regression <- size_regression

  model$build_expected_profile_and_sample_peak_heights <- function(genotypes){
    expected_profile <- log_normal_model_build_expected_profile(model, genotypes)

    x <- log_normal_model_sample_peak_heights(model,
                                              expected_profile,
                                              model$stutter_variability)

    x$LSAE <- LSAE[x$Marker]

    x
  }

  if (!missing(stutter_model)){
    model$stutter_model <- stutter_model
    model$stutter_variability <- stutter_variability
    model$parameters$k2 <- k2
  }

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

  degradation <- parameters$degradation

  x <- data.frame(
    Marker=character(),
    Allele=character(),
    Height=numeric(),
    Size=numeric(),
    Expected=numeric(),
    stringsAsFactors = FALSE)

  for (i_contributor in seq_len(number_of_contributors)){

    template_contributor <- parameters$template[i_contributor]

    g <- genotypes[[i_contributor]]

    for (i_row in seq_len(nrow(g))){

      locus <- g$Locus[i_row]
      ab <- c(g$Allele1[i_row], g$Allele2[i_row])

      lsae <- as.numeric(model$LSAE[locus])

      for (a in ab){
        size <- size_regression(locus, a)

        deg <- exp(-degradation[i_contributor] * (size - 100.))

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

  x$HeightAllele <- 10^(log10(x$ExpectedAllele) + rnorm(n = nrow(x),
                                                        mean = 0,
                                                        sd = sqrt(x$VarianceAllele)))

  i_stutter = 1

  # determine expected stutters based on _realised_ HeightAllele
  k2 <- parameters$k2

  for (i_stutter in seq_along(stutter_model$stutter_types)){
    stutter <- stutter_model$stutter_types[[i_stutter]]
    stutter_name <- stutter$name

    sr_column_name <- paste0("StutterRate", stutter_name)
    stutter_product_column_name <- paste0("StutterProduct", stutter_name)

    expected_column_name <- paste0("Expected", stutter_name)

    idx_parents <- !is.na(x[[stutter_product_column_name]])
    stutter_products <- x[[stutter_product_column_name]][idx_parents]

    idx_targets <- match(paste0(x$Marker[idx_parents], stutter_products),
                         paste0(x$Marker, x$Allele))

    x[[expected_column_name]][idx_targets] <- x$HeightAllele[idx_parents] *
                                    x$StutterRatioBackStutter[idx_parents]
  }

  # sample stutter heights
  stutter_variability <- model$stutter_variability

  stutter_types <- model$stutter_model$stutter_types

  for (stutter_name in names(stutter_types)){
    expected_column <- paste0("Expected", stutter_name)
    variance_column <- paste0("Variance", stutter_name)
    height_column <- paste0("Height", stutter_name)

    x[[height_column]] <- 0.

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

      x[[height_column]][idx_stutter] <- 10^(log10(x[[expected_column]][idx_stutter]) +
                                               rnorm(n = sum(idx_stutter),
                                                     mean = 0,
                                                     sd = sqrt(x[[variance_column]][idx_stutter])))
    }
    else{
      x[[variance_column]] <- stutter_k2 / (b / x[[expected_column]] + x[[expected_column]])

      x[[height_column]] <- 10^(log10(x[[expected_column]]) + rnorm(n = nrow(x),
                                                     mean = 0,
                                                     sd = sqrt(x[[variance_column]])))
    }
  }


  # add up stutters
  x$HeightStutter <- rowSums(x[paste0("Height", names(stutter_types))])

  x$Height <- x$HeightAllele + x$HeightStutter

  # fix expected stutter total
  x$ExpectedStutter <- rowSums(x[paste0("Expected", names(stutter_types))])

  x$Expected <- x$ExpectedStutter + x$ExpectedAllele


  x
}
