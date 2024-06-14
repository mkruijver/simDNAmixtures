.allele_tables_to_wide_references <- function(tabs){

  if (!is.list(tabs)) stop("Tabs need to be a list of allele tables")

  locus_names <- tabs[[1]]$Locus
  if (is.null(locus_names)) stop("A column named locus is needed")

  sample_names <- names(tabs)
  wide <- matrix(data = character(), nrow = length(sample_names),,
                 ncol = 2 * length(locus_names))
  rownames(wide) <- sample_names

  wide_rows <- lapply(tabs, function (y){
               data.frame(t(stats::setNames(
                 as.vector(rbind(y$Allele1, y$Allele2)),
                 paste0(rep(y$Locus, each=2),c("","(2)")))),
                 check.names = FALSE, stringsAsFactors = FALSE)
  })

  wide <- dplyr::bind_rows(wide_rows)



  rownames(wide) <- names(tabs)

  wide
}
