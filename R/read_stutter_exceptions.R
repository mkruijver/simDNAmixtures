#' @title Reads a stutter exceptions file with overrides for expected stutter ratios
#'
#' @param filename Character. Path to file.
#' @details Reads the file from disk and returns a named numeric vector with stutter ratio exceptions for a given locus and allele.
#' @examples
#' filename <- system.file("extdata","GlobalFiler_Stutter_Exceptions_3500.csv",package = "SimMixDNA")
#' exceptions <- read_stutter_exceptions(filename)
#' exceptions$TH01["9.3"]
#' @export
read_stutter_exceptions <- function(filename){

  exceptions_by_locus <- list()
  exceptions_table <- read.csv(filename)

  alleles <- exceptions_table$Allele

  loci <- names(exceptions_table)[-1]

  for (locus in loci){
    idx <- exceptions_table[[locus]] != 0

    exceptions_locus <- setNames(exceptions_table[[locus]][idx], alleles[idx])

    if (length(exceptions_locus) != 0){
      exceptions_by_locus[[locus]] <- exceptions_locus
    }
  }

  exceptions_by_locus
}