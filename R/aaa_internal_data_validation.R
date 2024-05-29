#' Validate Allele Frequencies
#'
#' This is an internal function to check the format and contents of the frequency list.
#'
#' @param freqs List of allele frequencies.
#' @param required_loci Character vector of required loci names.
#' @keywords internal
.validate_freqs <- function(freqs, required_loci = NULL){
  if (!is.list(freqs)){
    stop("freqs should be a list")
  }
  if (!all(sapply(freqs, is.numeric))){
    stop("freqs should be a list of numeric vectors")
  }

  if (!is.null(required_loci)){
    for (locus in loci){
      if (!(locus %in% names(freqs))){
        stop(paste0("freqs not available for locus "), locus)
      }
    }
  }
}
