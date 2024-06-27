#' @title Sample mixtures with known genotypes and random parameters according to priors
#'
#' @param n Integer. Number of samples.
#' @param genotypes List of contributor genotypes. See \link{sample_contributor_genotypes}.
#' @param sampling_parameters List. Passed to the sample_model function.
#' @param model_settings List. Passed to the sample_model function.
#' @param sample_model Function such as \link{sample_log_normal_model}.
#' @param results_directory (optionally) Character with path to directory where results are written to disk.
#' @param seed (optionally) Integer seed value that can be used to get reproducible runs. If results are written to disk, the 'Run details.txt' file will contain a seed that can be used for reproducing the result.
#' @param number_of_replicates Integer. Number of replicate simulations for each sample.
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
#' # define two known GlobalFiler genotypes
#' K1 <- data.frame(
#'   `Sample Name` = "K1",
#'   Locus = c("D3S1358", "vWA", "D16S539", "CSF1PO", "TPOX", "D8S1179",
#'             "D21S11", "D18S51", "D2S441", "D19S433", "TH01", "FGA",
#'             "D22S1045", "D5S818", "D13S317", "D7S820", "SE33",
#'             "D10S1248", "D1S1656", "D12S391", "D2S1338"),
#'   Allele1 = c("14", "15", "11", "11", "8",
#'               "11", "28", "14", "11", "14", "9.3",
#'               "21", "15", "9", "9", "10", "14",
#'               "14", "12", "21", "20"),
#'   Allele2 = c("16", "18", "13", "11", "11",
#'               "12", "30", "16", "12", "14", "9.3",
#'               "24", "16", "11", "12", "10", "18",
#'               "14", "15", "22", "24"), check.names = FALSE
#' )
#'
#' K2 <- data.frame(
#'   `Sample Name` = "K2",
#'   Locus = c("D3S1358", "vWA", "D16S539", "CSF1PO", "TPOX", "D8S1179",
#'             "D21S11", "D18S51", "D2S441", "D19S433", "TH01", "FGA",
#'             "D22S1045", "D5S818", "D13S317", "D7S820", "SE33",
#'             "D10S1248", "D1S1656", "D12S391", "D2S1338"),
#'   Allele1 = c("16", "16", "10", "13", "8",
#'               "11", "27", "17", "10", "13", "6",
#'               "23", "16", "11", "11", "9", "22.2",
#'               "14", "12", "22", "21"),
#'   Allele2 = c("16", "17", "12", "14", "8",
#'               "13", "30", "18", "11", "14", "6",
#'               "25", "17", "11", "11", "12", "28.2",
#'               "15", "14", "22", "23"), check.names = FALSE
#' )
#'
#' # first sample three replicates of a low-level profile of K1 only
#' gf <- gf_configuration()
#'
#' sampling_parameters <- list(min_template = 75., max_template = 75,
#'                             degradation_shape = 2.5, degradation_scale = 1e-3)
#'
#' single_source_results <- sample_mixtures_from_genotypes(n = 1,
#'                 genotypes = list(K1), sampling_parameters = sampling_parameters,
#'                 number_of_replicates = 3, sample_model = sample_log_normal_model,
#'                 model_settings = gf$log_normal_bwfw_settings)
#'
#' # now sample two mixtures of K1 and K2 with two replicates each
#'
#' mix_results <- sample_mixtures_from_genotypes(n = 2,
#'                 genotypes = list(K1, K2), sampling_parameters = sampling_parameters,
#'                 number_of_replicates = 3, sample_model = sample_log_normal_model,
#'                 model_settings = gf$log_normal_bwfw_settings)
#' @export
sample_mixtures_from_genotypes <- function(n,
                                           genotypes,
                                           sampling_parameters,
                                           model_settings,
                                           sample_model,
                                           results_directory,
                                           seed,
                                           number_of_replicates = 1L,
                                           tag = "simulation", silent = FALSE){

  # validate inputs
  .validate_integer(n, require_strictly_positive = TRUE)
  .validate_integer(number_of_replicates, require_strictly_positive = TRUE)

  # if seed is missing we obtain one here
  seed_validated <- .validate_or_generate_seed(seed)
  set.seed(seed_validated)

  number_of_contributors <- length(genotypes)

  write_to_disk <- FALSE
  if (!missing(results_directory)){
    write_to_disk <- TRUE

    results_dirs <- .init_results_directory(results_directory, tag,
                                            seed_validated, call = deparse(match.call()))
  }

  samples <- vector(mode = "list", length = n * number_of_replicates)

  sample_names <- character(n * number_of_replicates)

  # keep track of sample index even if replicates are used
  i_sample_out <- 1L
  for (i_sample in seq_len(n)){
    contributor_genotypes <- genotypes

    model <- sample_model(number_of_contributors = number_of_contributors,
                          sampling_parameters = sampling_parameters,
                          model_settings = model_settings)

    for (i_rep in seq_len(number_of_replicates)){
      replicate_label <- if (number_of_replicates > 1) paste0("_rep", i_rep) else ""

      sample_name <- paste0("sample", "_", sprintf("%04d", i_sample),
                            "_", model$sample_name_suffix, replicate_label)

      sample_names[i_sample_out] <- sample_name

      annotated_mixture <- sample_mixture_from_genotypes(contributor_genotypes, model, sample_name)

      mixture <- get_bare_mixture(annotated_mixture)

      samples[[i_sample_out]] <- list(sample_name = sample_name,
                                      contributor_genotypes = contributor_genotypes,
                                      model = model,
                                      annotated_mixture = annotated_mixture,
                                      mixture = mixture)

      if (write_to_disk){
        .write_mixture(samples[[i_sample_out]], results_dirs)

        ## knowns (wide table)
        write_knowns(contributor_genotypes, results_dirs$knowns_dir, sample_name)
      }

      i_sample_out <- i_sample_out + 1L
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
