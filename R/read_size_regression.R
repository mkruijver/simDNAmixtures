#' @title Reads a size regression file
#'
#' @param filename Character. Path to file.
#' @details Reads the file from disk and returns a function that provides the fragment length (bp) for a given locus and allele.
#' @examples
#' filename <- system.file("extdata","GlobalFiler_SizeRegression.csv",package = "SimMixDNA")
#' regression <- read_size_regression(filename)
#' regression("vWA", 12)
#' @export
read_size_regression <- function(filename){

  regression_df <- utils::read.csv(filename, colClasses = c("character", "numeric", "numeric"))
  regression_df_by_locus <- split(regression_df, regression_df$Locus)


  f <- function(locus, allele){

    regression_locus <- regression_df_by_locus[[locus]]

    if (is.null(regression_locus)){
      stop("No size regression available for locus ", locus)
    }

    intercept <- regression_locus$Intercept
    slope <- regression_locus$Slope

    intercept + slope * as.numeric(allele)
  }

  f
}
