get_parameter_summary <- function(samples){

  # collect all parameter names across samples and their max length
  parameter_names <- unique(c(lapply(samples, function(x)
                            names(x$model$parameters)), recursive=TRUE))

  get_max_length_by_parameter_name <- function(parameter_name, s){
    max(sapply(s, function(x) length(x$model$parameters[[parameter_name]])))
  }

  max_length_by_parameter_name <- sapply(parameter_names, function(parameter_name){
    get_max_length_by_parameter_name(parameter_name, samples)
  })

  # create a df for each sample
  dfs_by_sample <- list()
  for (i_sample in seq_along(samples)){
    # for this sample, collect dfs for each parameter
    parameters <- samples[[i_sample]]$model$parameters

    dfs_by_parameter_name <- list()

    for (parameter_name in parameter_names){
      values <- parameters[[parameter_name]]

      if (is.null(names(values))){
        # pad with NAs if necessary
        if (max_length_by_parameter_name[[parameter_name]] > length(values)){
          values <- c(values, rep(NA_real_,
                                  max_length_by_parameter_name[[parameter_name]] - length(values)))
        }

        if (max_length_by_parameter_name[[parameter_name]] > 1){
          dfs_by_parameter_name[[parameter_name]] <- stats::setNames(data.frame(t(values)),
                                                              paste0(parameter_name, seq_along(values)))
        }else{
          dfs_by_parameter_name[[parameter_name]] <- stats::setNames(data.frame(t(values)),
                                                              parameter_name)
        }
      }
      else{
        dfs_by_parameter_name[[parameter_name]] <- data.frame(t(values))
      }
    }

    names(dfs_by_parameter_name) <- NULL
    df_sample <- do.call(cbind, dfs_by_parameter_name)

    dfs_by_sample[[i_sample]] <- df_sample
  }

  parameter_summary <- data.frame(SampleName = names(samples),
                                  do.call(rbind, dfs_by_sample))

}
