# @param model_settings List. Possible parameters: \itemize{
#  \item locus_names. Character vector.
#  \item size_regression. Function, see \link{read_size_regression}.
#  }

.validate_drop_model_settings <- function(model_settings){

  locus_names <- model_settings$locus_names
  size_regression <- model_settings$size_regression

  .validate_character(locus_names)
}
