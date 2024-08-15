.sample_many_founders <- function(ped, freqs,
                                  unrelated_names = character(),
                                  number_of_replicates = 1L, loci = names(freqs),
                                  sex_locus_name = "AMEL",
                                  return_integer = FALSE){

  # if no pedigree is provided then we create a dummy of size 0
  if (missing(ped)){
    ped <- dummy_pedigree
  }

  number_of_pedigree_members <- length(ped$ID)
  number_of_extra_unrelateds <- length(unrelated_names)

  number_of_persons <- number_of_pedigree_members + number_of_extra_unrelateds
  number_of_loci <- length(loci)

  # set up an allele matrix for all persons for all replicates across all loci
  x <- matrix(data = NA_integer_,
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

    founder_allele_locus_cols_idx <- rep(c(2 * i_locus - 1, 2 * i_locus),
                                         times = number_of_founders * number_of_replicates)
    founder_allele_idx[, 2] <- founder_allele_locus_cols_idx

    is_sex_locus <- locus_name == sex_locus_name

    if (!is_sex_locus){
      # sample and put the alleles in the right places
      a_int <- sample.int(n = length(f_locus),
                          size = number_of_founder_alleles,
                          replace = TRUE, prob = f_locus)
      x[founder_allele_idx] <- a_int
    }
    else{
      .validate_sex_locus_freqs(a_locus, f_locus)

      pr_xy <- as.numeric(f_locus[AMEL_XY_packed])
      pr_xx <- as.numeric(f_locus[AMEL_XX_packed])

      # prepare a single chunk
      chunk <- integer(number_of_persons)
      ped_sexed_idx <- which(ped$SEX == 1 | ped$SEX == 2)
      chunk[ped_sexed_idx] <- ped$SEX[ped_sexed_idx]

      # repeat the single chunk
      chunk_repped <- rep(chunk, number_of_replicates)

      # fill in the blanks
      chunk_repped_blank_idx <- which(chunk_repped == 0L)
      chunk_repped[chunk_repped_blank_idx] <- sample.int(
        n = 2, size = length(chunk_repped_blank_idx),
        replace = TRUE, prob = c(pr_xy, pr_xx))

      x[, 2*i_locus - 1] <- chunk_repped
    }
  }

  if (return_integer) x else .many_genotypes_int_to_labels(x, freqs,
            sex_locus_name = sex_locus_name, loci = loci)
}
