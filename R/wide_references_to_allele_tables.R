.wide_references_to_allele_tables <- function(x){

  if (isTRUE(startsWith(tolower(names(x)[2]), "sample"))){
    rownames(x) <- x[[2]]
    return(.wide_references_to_allele_tables(x[-c(1, 2)]))
  }

  sample_names <- rownames(x)
  locus_names <- colnames(x)[seq(from = 1L, to = length(colnames(x)), by = 2L)]

  tabs <- vector(mode = "list", length = nrow(x))
  names(tabs) <- sample_names

  for (i in seq_len(nrow(x))){
    a1a2 <- matrix(unlist(x[i,]), ncol = 2, byrow = TRUE)
    colnames(a1a2) <- paste0("Allele", 1:2)

    tab <- data.frame("Sample Name" = sample_names[i],
               Locus = locus_names,
               a1a2, check.names = FALSE)

    tabs[[i]] <- tab
  }

  tabs
}
