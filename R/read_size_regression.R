#' @title Reads a stutter regression file
#'
#' @param filename Character. Path to file.
#' @details Reads the file from disk and returns a function that provides the fragment length (bp) for a given locus and allele.
#' @examples
#' filename <- system.file("extdata","GlobalFiler_SizeRegression.csv",package = "SimMixDNA")
#' regression <- read_size_regression(filename)
#' regression("vWA", 12)
#' @export
read_size_regression <- function(filename){

  regression_df <- read.csv(filename, colClasses = c("character", "numeric", "numeric"))
  regression_df_by_locus <- split(regression_df, regression_df$Locus)

  f <- function(locus, allele){
    intercept <- regression_df_by_locus[[locus]]$Intercept
    slope <- regression_df_by_locus[[locus]]$Slope

    intercept + slope * allele
  }

  f
}
