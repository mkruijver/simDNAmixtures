#' @title Sample mixture profile with provided genotypes
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

  if (!is.list(genotypes)){
    stop("genotypes is not a list of DataFrames")
  }
  if (!all(sapply(genotypes, is.data.frame))){
    stop("genotypes is not a list of DataFrames")
  }
  if (!is.character(sample_name)){
    stop("sample_name is not a character")
  }
  if (length(sample_name) != 1){
    stop("sample_name is not length 1")
  }

  profile <- model$build_expected_profile_and_sample_peak_heights(genotypes)

  profile <- reorder_profile(profile)

  profile$SampleName <- rep(sample_name, nrow(profile))

  # make SampleName the first column
  profile <- profile[c(length(profile), seq_len(length(profile) - 1))]

  # round size and height
  profile$Size <- round(profile$Size)
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
