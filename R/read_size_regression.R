#' @title Reads a size regression file
#'
#' @param filename Path to file (character).
#' @param exceptions Optionally a list providing sizes for alleles not covered by the regression. See examples for how this can be used to assign sizes to X and Y at the Amelogenin locus.
#' @param repeat_length_by_marker Optionally a named integer vector with repeat lengths by marker. If not provided, then a .3 allele will not convert to e.g. .75 for a tetranucleotide.
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
#'
#' # now add AMEL sizes
#' regression_with_AMEL <- read_size_regression(filename, exceptions = list(
#'                           AMEL = stats::setNames(c(98.5, 104.5),
#'                                                  nm = c("X", "Y"))))
#' # check that we can obtain size for X at AMEL
#' stopifnot(regression_with_AMEL("AMEL", "X") == 98.5)
#'
#' # pass in repeat_length_by_marker for more precise estimates
#' gf <- gf_configuration()
#'
#' regression_with_repeat_length <- read_size_regression(filename,
#'            repeat_length_by_marker = gf$repeat_length_by_marker)
#'
#' # obtain size for the 15.3 allele at the D1S1656 locus
#' stopifnot(regression_with_repeat_length("D1S1656", 15.3) ==
#'            121.628203912362 + 15.75 * 4.2170043570021)
#'
#' @export
read_size_regression <- function(filename, exceptions, repeat_length_by_marker){

  regression_df <- utils::read.csv(filename, colClasses = c("character", "numeric", "numeric"))
  regression_df_by_locus <- split(regression_df, regression_df$Locus)

  has_exceptions <- !missing(exceptions)
  has_repeat_length_by_marker <- !missing(repeat_length_by_marker)

  f <- function(locus, allele){

    # first check if there is an override (used for AMEL)
    if (has_exceptions){
      locus_exceptions <- exceptions[[locus]]

      if (!is.null(locus_exceptions)){
        if (allele %in% names(locus_exceptions)){
          size <- locus_exceptions[[allele]]
          return(size)
        }
      }
    }

    regression_locus <- regression_df_by_locus[[locus]]

    if (is.null(regression_locus)){
      stop("No size regression available for locus ", locus)
    }

    intercept <- regression_locus$Intercept
    slope <- regression_locus$Slope

    allele_numeric <- as.numeric(allele)

    # potentially convert the point repeats the right way if
    if (has_repeat_length_by_marker){
      repeat_length <- repeat_length_by_marker[[locus]]

      if (!is.null(repeat_length)){
        allele_numeric <- repeats_to_decimals(allele_numeric, repeat_length)
      }
    }

    intercept + slope * allele_numeric
  }

  f
}
