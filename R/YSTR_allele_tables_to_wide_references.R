.YSTR_allele_tables_to_wide_references <- function(tabs){

  if (!is.list(tabs)) stop("Tabs need to be a list of allele tables")

  locus_names <- tabs[[1]]$Locus
  if (is.null(locus_names)) stop("A column named locus is needed")

  sample_names <- names(tabs)

  wide_rows <- lapply(tabs, function (y){
    data.frame(t(stats::setNames(
      apply(y[startsWith(names(y), "Allele")], 1, function(x)
        {
        # collapse vector into single character
        paste0(x[!is.na(x)], collapse = ",")
        }), nm = locus_names)),
      check.names = FALSE, stringsAsFactors = FALSE)
  })

  wide <- dplyr::bind_rows(wide_rows)

  rownames(wide) <- names(tabs)

  wide
}
