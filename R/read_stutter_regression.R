#' @title Reads a stutter regression file
#'
#' @param filename Character. Path to file.
#' @param min_stutter_ratio Numeric.
#' @details Reads the file from disk and returns a function that provides the expected stutter ratio for a given locus and allele.
#' @examples
#' filename <- system.file("extdata","GlobalFiler_Stutter_3500.txt",package = "SimMixDNA")
#' regression <- read_stutter_regression(filename)
#' regression("vWA", 12)
#' @export
read_stutter_regression <- function(filename, min_stutter_ratio = 0.001){

  if (!is.numeric(min_stutter_ratio)){
    stop("min_stutter_ratio is not numeric")
  }
  if (length(min_stutter_ratio) != 1){
    stop("min_stutter_ratio does not have length 1")
  }

  regression_df <- read.csv(filename, colClasses = c("character", "numeric", "numeric"))
  regression_df_by_locus <- split(regression_df, regression_df$Locus)

  f <- function(locus, allele){

    regression_locus <- regression_df_by_locus[[locus]]

    if (is.null(regression_locus)){
      stop("No stutter regression available for locus ", locus)
    }

    intercept <- regression_locus$Intercept
    slope <- regression_locus$Slope

    max(min_stutter_ratio, intercept + slope * as.numeric(allele))
  }

  f
}
