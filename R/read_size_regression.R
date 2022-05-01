#' @title Reads a size regression file
#'
#' @param filename Path to file (character).
#' @details
#' Read a regression file from disk and returns a function that provides the fragment length (bp) for a given locus and allele.
#'
#' DNA profiles consist of the observed peaks (alleles or stutter products) at several loci as well as the peak heights and sizes. The size refers to the fragment length (bp). A linear relationship exists between the size of a peak and the size. When peaks are sampled in the \link{sample_mixture_from_genotypes} function, a size is assigned using a size regression. The \code{read_size_regression} function reads such a regression from disk.
#' @return A function that takes a locus name and allele as arguments and returns the size.
#' @examples
#' filename <- system.file("extdata",
#'                         "GlobalFiler_SizeRegression.csv",
#'                         package = "simDNAmixtures")
#'
#' regression <- read_size_regression(filename)
#'
#' # obtain size for the 12 allele at the vWA locus
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
