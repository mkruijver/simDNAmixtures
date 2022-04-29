#' @title Sample a genotype according to allele frequencies
#'
#' @param freqs Allele frquencies (see \link{read_allele_freqs})
#' @param loci Character vector of locus names (defaults to names attr. of \code{freqs})
#' @param label Sample name
#' @details A genotype is sampled randomly by drawing two alleles from allele frequencies.
#' @examples
#' # below we read an allele freqs and sample a genotype
#' filename <- system.file("extdata","FBI_extended_Cauc.csv",
#'                         package = "simDNAmixtures")
#' freqs <- read_allele_freqs(filename)
#' sample_genotype(freqs, loci = c("D3S1358", "vWA"))
#' @export
sample_genotype <- function(freqs, loci = names(freqs), label = "U"){

  if (!is.list(freqs)){
    stop("freqs should be a list")
  }
  if (!all(sapply(freqs, is.numeric))){
    stop("freqs should be a list of numeric vectors")
  }
  for (locus in loci){
    if (!(locus %in% names(freqs))){
      stop(paste0("freqs not available for locus "), locus)
    }
  }

  number_of_loci <- length(loci)

  epg <- data.frame("Sample Name" = rep(label, number_of_loci),
                    Locus = loci,
                    Allele1 = character(number_of_loci),
                    Allele2 = character(number_of_loci),
                    stringsAsFactors = FALSE, check.names = FALSE)

  for(i_locus in seq_along(loci)){
    locus <- loci[i_locus]

    f <- freqs[[locus]]

    tryCatch({
      ab <- names(f)[sort(sample.int(n = length(f),
                 size = 2, replace = TRUE, prob = f))]
    }, error = function(e){

      e$message <- paste0(e$message, " (at ", locus, ")")
      stop(e)
    })

    epg$Allele1[i_locus] <- ab[1]
    epg$Allele2[i_locus] <- ab[2]
  }

  epg
}
