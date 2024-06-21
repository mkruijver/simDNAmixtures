.sample_many_founders <- function(ped, freqs,
                                  unrelated_names = character(),
                                  number_of_replicates = 1L, loci = names(freqs),
                                  sex_locus_name = "AMEL"){

  # if no pedigree is provided then we create a dummy of size 0
  if (missing(ped)){
    ped <- dummy_pedigree
  }

  number_of_pedigree_members <- length(ped$ID)
  number_of_extra_unrelateds <- length(unrelated_names)

  number_of_persons <- number_of_pedigree_members + number_of_extra_unrelateds
  number_of_loci <- length(loci)

  # set up an allele matrix for all persons for all replicates across all loci
  x <- matrix(data = NA_character_,
              nrow = number_of_persons * number_of_replicates,
              ncol = 2 * number_of_loci)

  prefix <- if (number_of_replicates > 1) paste0("rep",
                                              rep(seq_len(number_of_replicates), each = number_of_persons), "_") else ""

  x_rownames <- paste0(prefix, rep(c(ped$ID, unrelated_names), number_of_replicates))
  rownames(x) <- x_rownames
  x_colnames <- paste0(rep(loci, each = 2), c("", "(2)"))
  colnames(x) <- x_colnames

  # sample founder alleles for all replicates
  idx_u_names <- if (number_of_extra_unrelateds==0) NULL else
     length(ped$FIDX) + seq_len(number_of_extra_unrelateds)
  ped_founder_row_idx <- c(which(ped$FIDX == 0L & ped$MIDX == 0L),
                           idx_u_names) # and extra unrelateds
  number_of_founders <- length(ped_founder_row_idx)

  replicate_row_offsets <- rep(seq(from = 0, to = number_of_replicates - 1),
                               each = 2 * number_of_founders) * number_of_persons

  founder_allele_rows_idx <- replicate_row_offsets +
    rep(rep(ped_founder_row_idx, each = 2), times = number_of_replicates)

  founder_allele_idx <- cbind(founder_allele_rows_idx, NA)

  number_of_founder_alleles <- 2 * number_of_founders * number_of_replicates

  for (i_locus in seq_along(loci)){
    locus_name <- loci[i_locus]

    f_locus <- freqs[[locus_name]]
    a_locus <- names(f_locus)

    if (length(f_locus) == 0){
      stop("Zero allele frequencies available for locus ", locus_name)
    }

    founder_allele_locus_cols_idx <- rep(c(2 * i_locus - 1, 2 * i_locus), times = number_of_founders * number_of_replicates)
    founder_allele_idx[, 2] <- founder_allele_locus_cols_idx

    is_sex_locus <- locus_name == sex_locus_name

    if (!is_sex_locus){
      # sample and put the alleles in the right places
      x[founder_allele_idx] <- a_locus[sample.int(n = length(f_locus),
                                                  size = number_of_founder_alleles,
                                                  replace = TRUE, prob = f_locus)]
    }
    else{
      .validate_sex_locus_freqs(a_locus, f_locus)

      # assign XX,XY to pedigree members and extra unrelated persons
      ped_sexed_idx <- which(ped$SEX == 1 | ped$SEX == 2)
      ped_unknown_sex_idx <- which(ped$SEX == 0)

      number_of_pedigree_members <- length(ped$SEX)

      ped_sexed_alleles <- sapply(ped$SEX[ped_sexed_idx], function(s)
                if (s==1) AMEL_XY_unpacked else AMEL_XX_unpacked)

      number_of_unsexed_persons <- number_of_extra_unrelateds + length(ped_unknown_sex_idx)

      # unpack alleles from provided genotypes
      a_locus_split <- lapply(a_locus, function(a) strsplit(a, split = ",")[[1]])
      extra_unrelateds_idx <- if (number_of_extra_unrelateds == 0) integer() else
        number_of_pedigree_members + seq_len(number_of_extra_unrelateds)

      for (i_rep in seq_len(number_of_replicates)){
        number_of_rows_offset <- (i_rep - 1) * number_of_persons

        # write XX,XY for this rep for the sexed ped members
        if (length(ped_sexed_alleles) > 0){
          x[cbind(rep(ped_sexed_idx, each = 2) +
                    number_of_rows_offset,
                  c(2 * i_locus - 1, 2 * i_locus))] <- ped_sexed_alleles
        }

        # sample XX,XY for unsexed ped members and unsexed unknown unrelated persons
        if (number_of_unsexed_persons > 0){
          genotypes_idx_sampled <- sample.int(n = length(f_locus),
                                              size = number_of_unsexed_persons,
                                              replace = TRUE, prob = f_locus)

          # grab the alleles for sampled genotypes
          unsexed_alleles <- unlist(a_locus_split[genotypes_idx_sampled])

          # put the XX,XY in the right places for these persons too
          x[cbind(rep(c(ped_unknown_sex_idx, extra_unrelateds_idx), each = 2) +
                  number_of_rows_offset,
                  c(2 * i_locus - 1, 2 * i_locus))] <- unsexed_alleles
        }
      }

    }
  }

  x
}
