

#  Possible parameters: \itemize{
#  \item locus_names. Character vector.
#  \item degradation_parameter_cap. Numeric.
#  \item c2_prior. Numeric of length two with shape and scale.
#  \item LSAE_variance_prior. Numeric of length one.
#  \item detection_threshold. Numeric vector (named) with Detection Thresholds. Defaults to 50 for each locus.
#  \item size_regression. Function, see \link{read_size_regression}.
#  \item stutter_model. Optionally a stutter_model object that gives expected stutter heights. See \link{global_stutter_model}.
#  \item stutter_variability. Optionally peak height variability parameters for stutters. Required when stutter_model is supplied.

validate_log_normal_model_settings <- function(model_settings, LSAE, k2,
                                               validate_LSAE = TRUE, validate_k2 = TRUE){

  locus_names <- model_settings$locus_names
  degration_parameter_cap <- model_settings$degradation_parameter_cap
  c2_prior <- model_settings$c2_prior
  LSAE_variance_prior <- model_settings$LSAE_variance_prior
  detection_threshold <- model_settings$detection_threshold
  size_regression <- model_settings$size_regression
  stutter_model <- model_settings$stutter_model
  stutter_variability <- model_settings$stutter_variability

  if (!is.character(locus_names)){
    stop("locus_names needs to be a character vector")
  }


  if (!all(locus_names %in% names(detection_threshold))){
    stop("all locus names need to be in names(detection_threshold)")
  }

  if (validate_LSAE){
    if (!is.numeric(LSAE)){
      stop("LSAE needs to be numeric")
    }

    if (!all(locus_names %in% names(LSAE))){
      stop("all locus names need to be in names(LSAE)")
    }
  }

  if (!is.numeric(detection_threshold)){
    stop("detection_threhold needs to be numeric")
  }

  if (!is.null(stutter_model)){
    if (is.null(stutter_variability)){
      stop("stutter_variability needs to be supplied when stutter_model is supplied")
    }

    if (validate_k2){
      if (missing(k2)){
        stop("k2 needs to be supplied when stutter_model is supplied")
      }

      expected_k2_names <- paste0("k2", names(stutter_model$stutter_types))

      if (!identical(expected_k2_names, names(k2))){
        stop("k2 does not have expected names: ", paste(expected_k2_names, collapse = ", "))
      }
    }

    for (stutter_name in names(stutter_variability)){
      if (!is.numeric(stutter_variability[[stutter_name]]$max_stutter_ratio)){
        stop("stutter_variability$", stutter_name,
             "$max_stutter_ratio is not numeric")
      }
      if (length(stutter_variability[[stutter_name]]$max_stutter_ratio) != 1){
        stop("stutter_variability$", stutter_name,
             "$max_stutter_ratio needs to have length 1")
      }
      if (stutter_variability[[stutter_name]]$max_stutter_ratio < 0){
        stop("stutter_variability$", stutter_name,
             "$max_stutter_ratio needs to be non-negative")
      }

    }
  }

}
