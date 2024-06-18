
# determines the minimum size of fragment lengths of all genotypes
# based on size regression
.get_deg_starts_at <- function(genotypes, size_regression){

  min_fragment_size <- Inf

  genotype <- genotypes[[1]]
  for (genotype in genotypes){
    allele_columns <- .get_allele_columns(genotype)

    for (i_locus in seq_len(nrow(allele_columns))){
      locus_name <- genotype$Locus[i_locus]

      for (i_allele in seq_len(ncol(allele_columns))){

        a <- allele_columns[[i_allele]][i_locus]

        if (!is.na(a)){
          size <- size_regression(locus_name, a)
          min_fragment_size <- min(min_fragment_size, size)
        }
      }
    }
  }

  min_fragment_size
}
