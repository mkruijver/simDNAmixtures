.many_genotypes_int_to_labels <- function(x, freqs, sex_locus_name,
                                          sort_alleles = TRUE,
                                          loci = names(freqs)){

  x_out <- matrix(data = NA_character_, nrow = nrow(x), ncol = ncol(x))
  dimnames(x_out) <- dimnames(x)

  for (i_locus in seq_along(loci)){
    locus_name <- loci[i_locus]

    f_locus <- freqs[[locus_name]]
    a_locus <- names(f_locus)

    is_sex_locus <- locus_name == sex_locus_name

    if (!is_sex_locus){
      a_int <- x[, 2*i_locus - 1]
      b_int <- x[, 2*i_locus]

      if (sort_alleles){
        a_character <- a_locus[pmin.int(a_int, b_int)]
        b_character <- a_locus[pmax.int(a_int, b_int)]
      }else{
        a_character <- a_locus[a_int]
        b_character <- a_locus[b_int]
      }

      x_out[, 2*i_locus - 1] <- a_character
      x_out[, 2*i_locus] <- b_character
    }
    else{
      x_out[, 2*i_locus - 1] <- "X"

      # note that male is 1 as is the ped format
      x_out[, 2*i_locus] <- c("Y", "X")[x[, 2*i_locus - 1]]
    }
  }

  x_out
}
