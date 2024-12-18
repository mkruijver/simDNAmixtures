#' Read STRmix Kit Settings
#'
#' Read STRmix kit settings from an XML file and associated stutter directory.
#'
#' @param filename Character vector specifying the path to the XML file containing the profiling kit settings.
#' @param stutters_dir Character vector specifying the directory path where stutter settings files are located.
#' @param include_y_loci Should Y loci be included in the list of locus names? Default is `FALSE`.
#'
#' @return A list containing the following components:
#' \itemize{
#'  \item locus_names. Character vector.
#'  \item degradation_parameter_cap. Numeric.
#'  \item c2_prior. Numeric of length two with shape and scale.
#'  \item LSAE_variance_prior. Numeric of length one.
#'  \item detection_threshold. Numeric vector (named) with Detection Thresholds. Defaults to 50 for each locus.
#'  \item size_regression. Function, see \link{read_size_regression}.
#'  \item stutter_model. A list representing the stutter model.
#'  \item stutter_variability. A list representing the stutter variability.
#' }
#'
#' @export
read_STRmix_kit_settings <- function(filename, stutters_dir,
                                     include_y_loci = FALSE){

  kit_xml <- filename |> xml2::read_xml() |> xml2::as_list()

  # check if this is a Y kit and warn if include_y_loci == FALSE
  kit_type <- kit_xml$profilingKit$kitType[[1]]
  if (isTRUE(grepl("YFiler|PowerPlex Y|PowerplexY|PPY23", kit_type, ignore.case = TRUE)) &&
      !include_y_loci){
    warning("include_y_loci = FALSE but kit type is ", kit_type)
  }

  # load stutters
  stutters <- .read_STRmix_kit_stutters(kit_xml, stutters_dir)

  size_regression_filename <- kit_xml$profilingKit$sizeRegressionFile[[1]]

  # read any size exceptions for AMEL
  profiling_kit_loci <- kit_xml$profilingKit$loci

  size_exceptions <- list()
  for (x in profiling_kit_loci){
    locus_name <- attr(x, "name")

    alleles <- x$alleles
    if (!is.null(alleles)){
      alleles_name <- sapply(alleles, function(a) attr(a, "name"))
      alleles_bp <- as.numeric(sapply(alleles, function(a) attr(a, "basePairs")))

      size_exceptions[[locus_name]] <- setNames(alleles_bp, alleles_name)
    }
  }

  # determine if there are any Y loci and quality loci
  quality_loci <- character()
  sex_loci <- character()
  y_loci <- character()

  for (x in profiling_kit_loci){
    locus_name <- attr(x, "name")

    if (.parse_STRmix_boolean(x$qualityMarker[[1]])){
      quality_loci <- c(quality_loci, locus_name)
    }
    if (.parse_STRmix_boolean(x$genderMarker[[1]])){
      sex_loci <- c(sex_loci, locus_name)
    }
    if (.parse_STRmix_boolean(x$yMarker[[1]])){
      y_loci <- c(y_loci, locus_name)
    }
  }


  repeat_length_by_marker <- .extract_repeat_length_by_locus_from_STRmix_kit(kit_xml)

  size_regression <- read_size_regression(
    file.path(stutters_dir, size_regression_filename),
    exceptions = size_exceptions,
    repeat_length_by_marker = repeat_length_by_marker)

  locus_names <- as.character(sapply(kit_xml$profilingKit$loci,
                                     function(x) attr(x, "name")))

  # ignore quality loci
  locus_names <- locus_names[!(locus_names %in% quality_loci)]
  if (!include_y_loci){
    locus_names <- locus_names[!(locus_names %in% y_loci)]
  }


  detection_threshold <- stats::setNames(sapply(kit_xml$profilingKit$kitSettings$detectionThresholds, function(x) as.numeric(x[[1]])),
                                  sapply(kit_xml$profilingKit$kitSettings$detectionThresholds, function(x) attr(x, "locus")))

  degradation_parameter_cap <- .parse_STRmix_double(kit_xml$profilingKit$kitSettings$degradationMax[[1]])

  c2_prior <- .parse_STRmix_double(kit_xml$profilingKit$kitSettings$allelicVariance[[1]])

  LSAE_variance_prior <- .parse_STRmix_double(kit_xml$profilingKit$kitSettings$locusAmpVariance[[1]])


  # create stutter model only if it is defined
  stutter_model <- NULL
  if (!is.null(stutters$stutter_model)){
    stutter_model <- allele_specific_stutter_model(stutters$stutter_model, size_regression)
  }

  # check that at least one locus is included
  if (length(locus_names) < 1){
    stop("No usable loci were found")
  }

  list(locus_names = locus_names,
       quality_loci = quality_loci,
       sex_loci = sex_loci,
       y_loci = y_loci,
       degradation_parameter_cap = degradation_parameter_cap,
       c2_prior = c2_prior,
       LSAE_variance_prior = LSAE_variance_prior,
       detection_threshold = detection_threshold,
       size_regression = size_regression,
       stutter_model = stutter_model,
       stutter_variability = stutters$stutter_variability)
}
