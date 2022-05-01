#' @title Sample genotypes for mixture contributors according to allele frequencies
#'
#' @param contributors Character vector with unique names of contributors. Valid names are "U1", "U2", ... for unrelated contributors or the names of pedigree members for related contributors.
#' @param freqs Allele frequencies (see \link{read_allele_freqs})
#' @param pedigree (optionally) [ped][pedtools::ped] object
#' @param loci Character vector of locus names (defaults to \code{names} attribute of \code{freqs})
#' @param return_non_contributors Logical. Should genotypes of non-contributing pedigree members also be returned?
#' @details For each founder, a genotype is sampled randomly by drawing two alleles from allele frequencies.
#' @examples
#'
#' # read allele frequencies
#' freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc.csv",
#'                            package = "simDNAmixtures"))
#'
#' # define a pedigree of siblings S1 and S2 (and their parents)
#' ped_sibs <- pedtools::nuclearPed(nch = 2,
#' father = "F", mother = "M",
#' children = c("S1", "S2"))
#'
#' # sample genotypes for a mixture of S1 + U1 + S2
#' # where U1 is an unrelated person
#'
#' sample_contributor_genotypes(contributors = c("S1","U1","S2"), freqs, ped_sibs)
#' @export
sample_contributor_genotypes <- function(contributors, freqs, pedigree,
                                         loci = names(freqs), return_non_contributors = FALSE){

  if (!is.character(contributors)){
    stop("contributors should be a character vector")
  }

  if (!is.logical(return_non_contributors)){
    stop("return_non_contributors should be a logical")
  }
  if (length(return_non_contributors) != 1){
    stop("return_non_contributors should be a logical of length 1")
  }

  if (!missing(pedigree)){
    if (!inherits(pedigree, "ped")){
      stop("pedigree should be of class ped")
    }

    ped_names <- pedigree$ID

    forbidden_names <- paste0("U", seq_len(10))

    forbidden_names_in_ped <- intersect(ped_names, forbidden_names)

    if (length(forbidden_names_in_ped) > 0){
      stop("Pedigree contains illegal name(s): ", paste(forbidden_names_in_ped, collapse = ", "))
    }

  }
  else{
    ped_names <- character()
  }

  # verify that names of other contributors are U1, U2, ..
  unr_names <- contributors[!(contributors %in% ped_names)]

  if (length(unr_names) > 0){
    expected_unr_names <- paste0("U", seq_along(unr_names))

    if (!setequal(expected_unr_names, unr_names)){
      stop("Expected unrelated contributor(s) named ", paste(expected_unr_names, collapse = ", "),
           " instead of ", paste(unr_names, collapse = ", ") )
    }
  }

  genotypes <- list()

  # sample related genotypes if any
  related_names <- contributors[contributors %in% ped_names]

  if (length(related_names) > 0){
    genotypes <- sample_pedigree_genotypes(pedigree = pedigree, freqs = freqs, loci = loci)
  }

  # sample unrelated genotypes
  for (unr_name in unr_names){
    genotypes[[unr_name]] <- sample_genotype(freqs = freqs, loci = loci, label = unr_name)
  }

  # reorder
  if (return_non_contributors){

    non_contributors_in_ped <- ped_names[!ped_names %in% contributors]

    genotypes <- genotypes[c(contributors, non_contributors_in_ped)]
  }
  else{
    genotypes <- genotypes[contributors]
  }

  genotypes
}
