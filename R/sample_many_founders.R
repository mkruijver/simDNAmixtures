.sample_many_founders <- function(ped, freqs,
                                  number_of_extra_unrelateds = 0L,
                                  number_of_replicates = 1L, loci = names(freqs)){

  number_of_persons <- length(ped$ID) + number_of_extra_unrelateds
  u_names <- if (number_of_extra_unrelateds==0) NULL else
    paste0("U", seq_len(number_of_extra_unrelateds))

  number_of_loci <- length(loci)

  # set up an allele matrix for all persons for all replicates across all loci
  x <- matrix(data = NA_character_,
              nrow = number_of_persons * number_of_replicates,
              ncol = 2 * number_of_loci)

  prefix <- if (number_of_replicates > 1) paste0("rep",
                                              rep(seq_len(number_of_replicates), each = number_of_persons), "
                                              _") else ""

  x_rownames <- paste0(prefix, rep(c(ped$ID, u_names), number_of_replicates))
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

    founder_allele_locus_cols_idx <- rep(c(2 * i_locus - 1, 2 * i_locus), times = number_of_founders * number_of_replicates)
    founder_allele_idx[, 2] <- founder_allele_locus_cols_idx

    x[founder_allele_idx] <- a_locus[sample.int(n = length(f_locus),
                                                size = number_of_founder_alleles,
                                                replace = TRUE, prob = f_locus)]
  }

  x
}
