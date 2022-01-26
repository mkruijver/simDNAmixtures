#' @title Read wide table (.txt) with Allele1, Allele2, ... columns as is
#'
#' @param filename Path to txt file.
#' @return Dataframe
#' @export
read_wide_table <- function(filename){
  x <- utils::read.csv(filename,sep = "\t",stringsAsFactors = FALSE,colClasses = "character",check.names = FALSE)
  stopifnot(all(sapply(x,class)=="character"))

  # determine columns per locus
  allele_cols <- startsWith(colnames(x),"Allele")
  size_cols <- startsWith(colnames(x),"Size")
  height_cols <- startsWith(colnames(x),"Height")

  can_import_sizeCols <- sum(allele_cols)==sum(size_cols)
  if (!can_import_sizeCols) warning("Could not read allele sizes")
  can_import_heightCols <- sum(allele_cols)==sum(height_cols)
  if (!can_import_heightCols) warning("Could not read peak heights")

  # convert sizes/heights to numeric
  cols <- colnames(x)[size_cols | height_cols]

  for (col in cols){
    x[[col]] <- as.numeric(x[[col]])
  }

  x
}
