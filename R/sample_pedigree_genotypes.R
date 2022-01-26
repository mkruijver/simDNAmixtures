#' @title Sample genotypes for pedigree according to allele frequencies
#'
#' @param pedigree Ped object
#' @param freqs Allele frequencies (see \link{read_allele_freqs})
#' @param loci Character vector of locus names (defaults to names attr. of \code{freqs})
#' @details For each founder, a genotype is sampled randomly by drawing two alleles from allele frequencies.
#' @examples
#' freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc.csv",package = "SimMixDNA"))
#'
#' ped_sibs <- pedtools::nuclearPed(nch = 2,
#' father = "F", mother = "M",
#' children = c("S1", "S2"))
#' sibs_genotypes <- SimMixDNA::sample_pedigree_genotypes(ped = ped_sibs,
#' freqs = freqs, loci = get_GlobalFiler_3500_data()$autosomal_markers)
#' @export
sample_pedigree_genotypes <- function(pedigree, freqs, loci = names(freqs)){

  if (!inherits(pedigree, "ped")){
    stop("pedigree should be of class ped")
  }
  if (!is.list(freqs)){
    stop("freqs should be a list")
  }
  if (!all(sapply(freqs, is.numeric))){
    stop("freqs should be a list of numeric vectors")
  }
  for (locus in loci){
    if (!(locus %in% names(freqs))){
      stop(paste0("freqs not available for locus "), locus)
    }
  }

  # sample founders
  profiles_by_id <- list()

  ped_row_is_founder <- pedigree$FIDX==0 & pedigree$MIDX==0

  founder_ids <- pedigree$ID[ped_row_is_founder]
  for (founder_id in founder_ids){
    profiles_by_id[[founder_id]] <- sample_genotype(freqs, loci, label = founder_id)
  }

  ped_row_is_sampled <- ped_row_is_founder

  # keep sampling non-founders whose both parents are sampled
  # until all non-founders have a genotype

  # sample non-founders
  while(!all(ped_row_is_sampled)){
    for (i_row in which(!ped_row_is_sampled)){
      person_id <- pedigree$ID[i_row]
      father_id <- pedigree$ID[pedigree$FIDX[i_row]]
      mother_id <- pedigree$ID[pedigree$MIDX[i_row]]

      if (!(is.null(profiles_by_id[[father_id]]) ||
          is.null(profiles_by_id[[mother_id]]))){

        profiles_by_id[[person_id]] <-
          sample_offspring(father = profiles_by_id[[father_id]],
                           mother = profiles_by_id[[mother_id]],
                           label = person_id)

        ped_row_is_sampled[i_row] <- TRUE
      }
    }
  }

  profiles_by_id
  # profiles <- do.call(rbind, profiles_by_id[pedigree$ID])
  # rownames(profiles) <- NULL
  #
  # profiles
}
