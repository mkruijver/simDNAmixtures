unpack_haplotypes <- function(haplotypes){
  strsplit(haplotypes,split = ",")
}

.wide_YSTR_references_to_allele_tables <- function(x){

  sample_names <- rownames(x)
  locus_names <- names(x)

  tabs <- vector(mode = "list", length = nrow(x))
  names(tabs) <- sample_names

  # split such that e.g. "12,13" becomes c("12", "13")
  x_unpacked_by_locus <- sapply(x, unpack_haplotypes)
  x_unpacked_by_locus_length <- sapply(x_unpacked_by_locus, length)
  max_number_of_alleles <- max(x_unpacked_by_locus_length)

  # process profile by profile (i.e. row by row)
  for (i_row in seq_len(nrow(x))){

    x_matrix <- matrix(data = character(),
                       nrow = length(locus_names),
                       ncol = max_number_of_alleles)

    ref_unpacked <- x_unpacked_by_locus[i_row,]

    for (i_locus in seq_along(locus_names)){
      a <- ref_unpacked[[i_locus]]
      x_matrix[i_locus, ][1:length(a)] <- a
    }

    colnames(x_matrix) <- paste0("Allele", seq_len(max_number_of_alleles))

    tab <- data.frame("Sample Name" = sample_names[i_row],
                      Locus = locus_names,
                      x_matrix, check.names = FALSE)

    tabs[[i_row]] <- tab
  }

  tabs
}
