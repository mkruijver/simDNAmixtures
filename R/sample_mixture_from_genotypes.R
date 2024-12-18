#' @title Sample mixture profile with provided genotypes
#'
#' @param genotypes List of contributor genotypes. See \link{sample_contributor_genotypes}.
#' @param model pg_model object.
#' @param sample_name Character. Defaults to "mixture".
#' @details A mixture profile is sampled according to the provided \code{pg_model} (see \link{gamma_model}, \link{log_normal_model} and genotypes (see \link{sample_contributor_genotypes}).
#' @return DataFrame with at least SMASH columns (see \link{SMASH_to_wide_table}). Depending on the chosen \code{pg_model} (e.g. \link{gamma_model} or \link{log_normal_model}), other columns with further details about the simulation are returned as well.
#' @seealso \link{sample_mixtures} for a function that samples many mixtures in one go.
#' @examples
#' # read allele frequencies and kit data
#' freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc_022024.csv",
#'                            package = "simDNAmixtures"))
#' gf <- gf_configuration()
#'
#' # define a pedigree of siblings S1 and S2 (and their parents)
#' ped_sibs <- pedtools::nuclearPed(children = c("S1", "S2"))
#'
#' # sample genotypes for a mixture of S1 + U1 + S2
#' # where U1 is an unrelated person
#' genotypes <- sample_contributor_genotypes(contributors = c("S1","U1","S2"),
#' freqs, pedigree = ped_sibs, loci = gf$autosomal_markers)
#'
#' # define a gamma model for peak heights
#' gamma_model <- gamma_model(mixture_proportions = c(0.5, 0.3, 0.2), mu = 1000.,
#'                     cv = 0.1, model_settings = gf$gamma_settings_no_stutter)
#'
#' # sample mixture from genotypes
#' mix <- sample_mixture_from_genotypes(genotypes, gamma_model)
#' @export
sample_mixture_from_genotypes <- function(genotypes, model, sample_name = "mixture"){

  if (!is.list(genotypes)){
    stop("genotypes is not a list of DataFrames")
  }
  if (!all(sapply(genotypes, is.data.frame))){
    stop("genotypes is not a list of DataFrames")
  }
  .validate_character(sample_name, required_length = 1L)
  .validate_character(model$locus_names, required_length_min = 1L)

  profile <- model$build_expected_profile_and_sample_peak_heights(genotypes)

  profile <- reorder_profile(profile)

  profile$SampleName <- rep(sample_name, nrow(profile))

  # make SampleName the first column
  profile <- profile[c(length(profile), seq_len(length(profile) - 1))]

  # round size and height
  profile$Size <- round(profile$Size, digits = 2)
  profile$Height <- round(profile$Height)

  profile
}

reorder_profile <- function(x){

  marker_order <- unique(x$Marker)

  # first order by size
  x1 <- x[order(x$Size),]

  # then by locus
  x2 <- x1[order(match(x1$Marker, table = marker_order)),]

  rownames(x2) <- NULL
  x2
}
