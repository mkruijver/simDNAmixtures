#' @title Read allele frequencies in STRmix a.k.a FSIgen format (.csv)
#'
#' @param filename Path to csv file.
#' @details Reads allele frequencies from a .csv file. The file should be in FSIgen format, i.e.
#'          comma separated with the first column specifying the allele labels and
#'          one column per locus. The last row should be the number of observations.
#'          No error checking is done since the file format is only loosely defined,
#'          e.g. we do not restrict the first column name or the last row name.
#' @return list
#' @examples
#' # below we read an allele freqs file that comes with the package
#' filename <- system.file("extdata","FBI_extended_Cauc.csv",package = "SimMixDNA")
#' freqs <- read_allele_freqs(filename)
#' freqs # the output is just a list with an N attribute
#' @importFrom utils read.csv
#' @export
read_allele_freqs <- function(filename){
  raw <- readLines(filename)

  dfWithoutN <- read.csv(file = filename,header = TRUE,nrows = length(raw)-2, check.names=FALSE)
  dfWithN <- read.csv(file = filename,header = TRUE,nrows = length(raw)-1, check.names=FALSE)

  returnList <- list()
  alleles <- dfWithoutN[[1]]

  indicesLoci <- seq(from=2,to=ncol(dfWithoutN),by=1)
  N <- setNames(numeric(length(indicesLoci)),nm = names(dfWithoutN[-1]))
  for(iLocus in indicesLoci){
    f0 <- dfWithoutN[[iLocus]]
    f0[is.na(f0)] <- 0.
    f <- f0

    returnList[[names(dfWithoutN)[iLocus]]]  <- setNames(f,nm =  alleles)
    N[names(dfWithoutN)[iLocus]] <- dfWithN[[iLocus]][length(dfWithN[[iLocus]])]
  }
  attr(returnList,"N") <- N

  returnList
}
