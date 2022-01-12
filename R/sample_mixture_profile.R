#' @title Sample mixture profile
#'
#' @param genotypes List of contributor genotypes. See \link{sample_contributor_genotypes}.
#' @param model pg_model object.
#' @param sample_name Character. Defaults to "mixture".
#' @details A mixture profile is sampled according to the provided pg_model.
#' @examples
#' # read allele frequencies
#' freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc.csv",package = "SimMixDNA"))
#'
#' # define a pedigree of siblings S1 and S2 (and their parents)
#' ped_sibs <- pedtools::nuclearPed(nch = 2,
#' father = "F", mother = "M",
#' children = c("S1", "S2"))
#'
#' # sample genotypes for a mixture of S1 + U1 + S2
#' # where U1 is an unrelated perso
#' genotypes <- sample_contributor_genotypes(contributors = c("S1","U1","S2"), ped_sibs, freqs)#'
#' @export
sample_mixture_profile <- function(genotypes, model, sample_name = "mixture"){

  profile <- model$build_expected_profile(genotypes)
  profile <- model$sample_peak_heights(profile)

  profile <- reorder_profile(profile)

  profile$SampleName <- rep(sample_name, nrow(profile))

  # make SampleName the first column
  profile <- profile[c(length(profile), seq_len(length(profile) - 1))]

  profile
}

reorder_profile <- function(x){

  locus_order <- unique(x$Locus)

  # first order by size
  x1 <- x[order(x$Size),]

  # then by locus
  x2 <- x1[order(match(x1$Locus, table = locus_order)),]

  rownames(x2) <- NULL
  x2
}
