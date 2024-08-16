#' @title Sample (many) genotypes for pedigree members according to allele frequencies by allele dropping and possibly taking linkage into account by simulating recombination.
#'
#' @param pedigree [ped][pedtools::ped] object
#' @param freqs Allele frequencies (see \link{read_allele_freqs})
#' @param loci Character vector of locus names (defaults to \code{names} attribute of \code{freqs})
#' @param unrelated_names Character vector with names of any additional unrelated persons. Defaults to length zero.
#' @param linkage_map A linkage map specifying the recombination fractions between loci. If NULL, loci are assumed to be independent.
#' @param number_of_replicates An integer specifying the number of replicate genotype samples to generate. Defaults to 1.
#' @param sex_locus_name Character vector, defaults to "AMEL"
#'
#' @seealso \link{read_allele_freqs}
#' @export
sample_many_pedigree_genotypes <- function(pedigree, freqs, loci = names(freqs),
                                           unrelated_names = character(),
                                           linkage_map,
                                           number_of_replicates = 1L,
                                           sex_locus_name = "AMEL"){

  if (missing(pedigree)){
    pedigree <- dummy_pedigree
  }

  .validate_pedigree(pedigree, disallow_U_names = TRUE)
  .validate_freqs(freqs, loci)

  locus_idx_by_name <- stats::setNames(seq_along(loci), loci)

  if (missing(linkage_map)) linkage_map <- NULL
  if (!is.null(linkage_map)) .validate_linkage_map(linkage_map)

  # add missing loci to linkage map (but not AMEL)
  missing_loci <- loci[!loci %in% c(linkage_map$locus, sex_locus_name)]
  linkage_map_used <- rbind(linkage_map[linkage_map$locus %in% loci,],
                            data.frame(chromosome = rep("missing", length(missing_loci)),
                                       locus = missing_loci,
                                       position = rep(NA, length(missing_loci))))

  # sample founders
  x <- .sample_many_founders(pedigree, number_of_replicates = number_of_replicates,
                             unrelated_names = unrelated_names,
                             freqs = freqs, loci = loci,
                             sex_locus_name = sex_locus_name,
                             return_integer = TRUE)

  # if there are no non-founders, then we are done here!

  # prepare indices for dropping alleles
  ped_row_is_founder <- pedigree$FIDX == 0 & pedigree$MIDX == 0
  ped_row_is_non_founder <- !ped_row_is_founder

  if (!any(ped_row_is_non_founder)){
    return(.many_genotypes_int_to_labels(x, freqs = freqs,
            sex_locus_name = sex_locus_name, loci = loci))
  }

  ped_non_founder_row_idx <- which(ped_row_is_non_founder)
  ped_non_founder_fidx <- pedigree$FIDX[ped_non_founder_row_idx]
  ped_non_founder_midx <- pedigree$MIDX[ped_non_founder_row_idx]

  transmissions <- data.frame(non_founder_row = rep(ped_non_founder_row_idx, each = 2),
                              fidx = rep(ped_non_founder_fidx, each = 2),
                              midx = rep(ped_non_founder_midx, each = 2),
                              allele = if (length(ped_non_founder_row_idx)>0) c(1, 2) else integer())

  number_of_persons <- length(pedigree$ID) + length(unrelated_names)

  replicate_row_offsets <- rep(seq(from = 0, to = number_of_replicates - 1),
                               each = nrow(transmissions)) * number_of_persons

  # for every locus we will adjust column 2 here
  to_idx <- cbind(rep(transmissions$non_founder_row, times = number_of_replicates) + replicate_row_offsets,
                  NA)
  # for every locus we adjust column 2 here based on the transmission_vectors and the locus idx
  from_idx <- cbind(rep(ifelse(transmissions$allele == 1L,
                               transmissions$midx,
                               transmissions$fidx),
                        times = number_of_replicates) +
                      replicate_row_offsets,
                    NA)

  # split the linkage map by chromosome
  linkage_map_by_chromosome <- split(linkage_map_used, linkage_map_used$chromosome)
  chromosomes <- naturalsort::naturalsort(names(linkage_map_by_chromosome))

  chromosome = chromosomes[1]
  # start sampling data by chromosome
  for (chromosome in chromosomes){

    # prepare linkage map for this chromosome
    linkage_map_chromosome <- linkage_map_by_chromosome[[chromosome]]
    linkage_map_chromosome$recombination_rate <- NA_real_

    for (i_row in seq_len(nrow(linkage_map_chromosome))[-1]){
      delta_position <- linkage_map_chromosome$position[i_row] - linkage_map_chromosome$position[i_row - 1]

      linkage_map_chromosome$recombination_rate[i_row] <- if (is.na(delta_position))
        0.5 else pedprobr::haldane(cM = delta_position)
    }

    # sample locus-by-locus
    chromosome_i_locus <- 1L
    for (chromosome_i_locus in seq_len(nrow(linkage_map_chromosome))){
      # sample a starting transmission vector for each replicate
      if (chromosome_i_locus == 1){
        transmission_vectors <- matrix(sample.int(n = 2,
                                                  size = nrow(transmissions) * number_of_replicates,
                                                  replace = TRUE),
                                       nrow = nrow(transmissions))
      }else{
        recombination_rate <- linkage_map_chromosome$recombination_rate[chromosome_i_locus]

        swapped <- c(2,1)[transmission_vectors]
        swap <- sample(c(TRUE, FALSE),
                       size = length(transmission_vectors), replace = TRUE,
                       prob = c(recombination_rate, 1.0 - recombination_rate))

        transmission_vectors <- matrix(ifelse(swap, yes = swapped, no = transmission_vectors),
                                       nrow = nrow(transmissions))
      }

      # determine the index of the loci in the output
      locus_idx <- as.integer(locus_idx_by_name[linkage_map_chromosome$locus[chromosome_i_locus]])
      # drop alleles down the pedigree for this locus
      from_idx[, 2] <- as.vector(transmission_vectors) + (2 * (locus_idx - 1))
      to_idx[, 2] <- rep(transmissions$allele, times = number_of_replicates) + (2 * (locus_idx - 1))

      for (i in seq_len(nrow(to_idx))){
        x[to_idx[i, , drop = FALSE]] <- x[from_idx[i, , drop=FALSE]]
      }
    }
  }

  return(.many_genotypes_int_to_labels(x, freqs = freqs,
                                       sex_locus_name = sex_locus_name, loci = loci))
}
