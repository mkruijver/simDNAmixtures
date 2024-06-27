#' @title Sample mixtures with known genotypes and fixed parameters
#'
#' @param genotypes List of contributor genotypes. See \link{sample_contributor_genotypes}.
#' @param parameter_summary DataFrame with parameters for each sample.
#' @param model_settings List. Passed to the sample_model function.
#' @param results_directory (optionally) Character with path to directory where results are written to disk.
#' @param seed (optionally) Integer seed value that can be used to get reproducible runs. If results are written to disk, the 'Run details.txt' file will contain a seed that can be used for reproducing the result.
#' @param tag Character. Used for sub directory name when results_directory is provided.
#' @param silent Logical. If TRUE, then no message will be printed about where the output (if any) was written to disk.
#' @return If \code{results_directory} is provided, this function has the side effect of writing results to disk.
#'
#' Return value is a list with simulation results:\itemize{
#' \item \code{call} matched call
#' \item \code{smash} DataFrame with all samples in SMASH format (see \link{SMASH_to_wide_table})
#' \item \code{samples} Detailed results for each sample
#' \item \code{parameter_summary} DataFrame with parameters for each sample
#' }
#' @seealso [sample_mixtures]
#' @examples
#' # simulate autosomal samples with fixed parameters (from csv) and refs (from csv)
#' parameter_summary <- utils::read.csv(system.file(
#'"extdata","Example_2p_Parameter_Summary.csv",
#'package = "simDNAmixtures"))
#'
#'gf <- gf_configuration()
#'
#'filename_refs <- system.file("extdata","Example_References_DB.csv",
#'                             package = "simDNAmixtures")
#'db_refs <- utils::read.csv(filename_refs, check.names = FALSE)
#'
#'genotypes <- simDNAmixtures:::.wide_references_to_allele_tables(db_refs)
#'
#'samples <- sample_mixtures_fixed_parameters(genotypes = genotypes,
#'                                            parameter_summary = parameter_summary,
#'                                            model_settings = gf$log_normal_bwfw_settings,
#'                                            seed = 1)
#' @export
sample_mixtures_fixed_parameters <- function(genotypes,
                                           parameter_summary,
                                           model_settings,
                                           results_directory,
                                           seed,
                                           tag = "simulation", silent = FALSE){

  # if seed is missing we obtain one here
  seed_validated <- .validate_or_generate_seed(seed)
  set.seed(seed_validated)

  # number_of_contributors <- length(genotypes)

  write_to_disk <- FALSE
  if (!missing(results_directory)){
    write_to_disk <- TRUE

    results_dirs <- .init_results_directory(results_directory, tag,
                                            seed_validated, call = deparse(match.call()))
  }

  number_of_samples <- nrow(parameter_summary)
  samples <- vector(mode = "list", length = number_of_samples)

  sample_names <- parameter_summary$SampleName

  # keep track of sample index even if replicates are used
  i_sample_out <- 1L

  contributors_columns <- which(startsWith(names(parameter_summary), "contributors"))

  get_parameter <- function(parameter_summary, i_sample, parameter_name){
    i_columns <- which(startsWith(names(parameter_summary), parameter_name))

    unlist(stats::na.omit(parameter_summary[i_sample, i_columns]))
  }

  for (i_sample in seq_len(number_of_samples)){

    contributor_names <- stats::na.omit(as.character(
      parameter_summary[i_sample, contributors_columns]))

    contributor_genotypes <- genotypes[contributor_names]

    model_type <- parameter_summary$model[i_sample]

    if (model_type == "log_normal_model"){
      if (!is.null(model_settings$stutter_model)){
        model <- log_normal_model(
          template = get_parameter(parameter_summary, i_sample, "template"),
          degradation = get_parameter(parameter_summary, i_sample, "template"),
          LSAE = unlist(parameter_summary[i_sample,model_settings$locus_names]),
          c2 = get_parameter(parameter_summary, i_sample, "c2"),
          k2 = get_parameter(parameter_summary, i_sample, "k2"),
          model_settings = model_settings)
      }else{
        model <- log_normal_model(
          template = get_parameter(parameter_summary, i_sample, "template"),
          degradation = get_parameter(parameter_summary, i_sample, "template"),
          LSAE = unlist(parameter_summary[i_sample,model_settings$locus_names]),
          c2 = get_parameter(parameter_summary, i_sample, "c2"),
          model_settings = model_settings)
      }

    }
    else if (model_type == "gamma_model"){
      stop("gamma_model is not implemented yet")
    }
    else{
      stop("unknown model type:", model_type)
    }

    sample_name <- sample_names[i_sample]

    annotated_mixture <- sample_mixture_from_genotypes(contributor_genotypes,
                                                       model, sample_name)

    mixture <- get_bare_mixture(annotated_mixture)

    samples[[i_sample]] <- list(sample_name = sample_name,
                                    contributor_genotypes = contributor_genotypes,
                                    model = model,
                                    annotated_mixture = annotated_mixture,
                                    mixture = mixture)

    if (write_to_disk){
      .write_mixture(samples[[i_sample]], results_dirs)

      ## knowns (wide table)
      write_knowns(contributor_genotypes, results_dirs$knowns_dir, sample_name)
    }
  }

  names(samples) <- sample_names

  # prepare additional outputs
  summaries <- .prepare_summaries(samples)

  if (write_to_disk){
    .write_summaries(summaries, tag, results_dirs)

    write(c(paste0("Simulation finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
          file = results_dirs$run_details_file, append = TRUE)

    if (!silent) cat("Finished sampling. Output written to", results_dirs$sub_dir, "\n")
  }

  list(call = match.call(),
       samples = samples,
       smash = summaries$smash,
       parameter_summary = summaries$parameter_summary)
}
