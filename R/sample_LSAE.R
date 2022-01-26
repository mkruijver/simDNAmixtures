#' @title Sample Locus Specific Amplification Efficiency (LSAE) according to prior
#'
#' @param LSAE_variance Numeric. See \link{get_GlobalFiler_3500_data} for an example.
#' @param locus_names Character vector.
#' @details In the Bright et al. log normal model, the expected peak height
#' includes a multiplicative factor for the locus (marker). These factors are
#' called the LSAEs (Locus Specific Amplification Efficiencies). In the model,
#' the prior for the log10 of LSAE is normal with mean 0. The variance
#' can be specified.
#' @examples
#'  gf <- get_GlobalFiler_3500_data()
#'  lsae <- sample_LSAE(gf$log_normal_settings$LSAE_variance_prior, gf$autosomal_markers)
#' @export
sample_LSAE <- function(LSAE_variance, locus_names){

  if (!is.numeric(LSAE_variance)){
    stop("LSAE_variance needs to be numeric")
  }

  if (length(LSAE_variance) != 1){
    stop("LSAE variance needs to have length 1")
  }

  if (LSAE_variance < 0){
    stop("LSAE_variance needs to be non-negative")
  }

  if (!is.character(locus_names)){
    stop("locus_names needs to be a character vector")
  }

  LSAE <- 10 ^ stats::rnorm(n = length(locus_names),
                   mean = 0,
                   sd = sqrt(LSAE_variance))

  stats::setNames(LSAE, locus_names)
}
