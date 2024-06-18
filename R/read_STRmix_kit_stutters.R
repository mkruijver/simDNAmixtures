.read_STRmix_kit_stutters <- function(kit_xml, stutters_dir){
  if (!is.list(kit_xml$profilingKit)){
    stop("kit_xml needs to be a list with an element named profilingKit")
  }

  # load general settings
  locus_names <- as.character(sapply(kit_xml$profilingKit$loci, function(x) attr(x,"name")))

  repeat_length_by_locus <- setNames(
    sapply(kit_xml$profilingKit$loci, function(x) as.numeric(x$repeatLength[[1]])),
    sapply(kit_xml$profilingKit$loci, function(x) attr(x, "name")))

  stutters <- list()
  stutter_variability <- list()

  stutters_settings <- kit_xml$profilingKit$kitSettings$stutters
  for(stutter_settings in stutters_settings){
    is_enabled <- .parse_STRmix_boolean(stutter_settings$enabled[[1]])

    if (!is_enabled){
      next
    }

    stutter_name <- attr(stutter_settings, "name")

    position_relative_to_parent <- .parse_STRmix_double(stutter_settings$positionRelativeToParent[[1]])


    regression_file <- stutter_settings$regressionFile[[1]]
    exceptions_file <- if (length(stutter_settings$exceptionsFile) == 0) NULL else
      stutter_settings$exceptionsFile[[1]]

    applicable_loci <- as.character(unlist(stutter_settings$applicableLoci))

    applies_to_all_loci <- all(locus_names %in% applicable_loci)

    stutter_regression <- read_stutter_regression(file.path(stutters_dir, regression_file))
    stutter_exceptions <- if (is.null(exceptions_file)) NULL else read_stutter_exceptions(file.path(stutters_dir, exceptions_file))

    # build the stutter model
    if (is.null(stutter_exceptions)) {
      stutters[[stutter_name]] <- stutter_type(name = stutter_name,
                                               delta = position_relative_to_parent,
                   repeat_length_by_marker = repeat_length_by_locus,
                   stutter_regression = stutter_regression,
                   applies_to_all_loci = applies_to_all_loci,
                   applies_to_loci = applicable_loci)
    } else {
      stutters[[stutter_name]] <- stutter_type(name = stutter_name,
                                               delta = position_relative_to_parent,
                   repeat_length_by_marker = repeat_length_by_locus,
                   stutter_regression = stutter_regression,
                   stutter_exceptions = stutter_exceptions,
                   applies_to_all_loci = applies_to_all_loci,
                   applies_to_loci = applicable_loci)
    }

    # and write out the stutter variability parameters
    inversely_proportional_to_parent <- .parse_STRmix_boolean(
      stutter_settings$inverselyProportionalToParent[[1]])

    stutter_max <- .parse_STRmix_double(stutter_settings$stutterMax[[1]])

    stutter_variance_prior <- .parse_STRmix_double(stutter_settings$stutterVariance[[1]])

    stutter_variability[[stutter_name]] <- list(k2_prior = stutter_variance_prior,
                                                inversely_proportional_to_parent = inversely_proportional_to_parent,
                                                max_stutter_ratio = stutter_max)

  }

  list(stutter_model = stutters,
       stutter_variability = stutter_variability)
}
