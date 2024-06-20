#' @title Sample mixtures with random genotypes and random parameters according to priors
#'
#' @param n Integer. Number of samples.
#' @param contributors Character vector with unique names of contributors. Valid names are "U1", "U2", ... for unrelated contributors or the names of pedigree members for related contributors.
#' @param freqs Allele frequencies (see \link{read_allele_freqs})
#' @param linkage_map (optional) A linkage map specifying the recombination fractions between loci. If missing, loci are assumed to be independent. See also \link{sample_many_pedigree_genotypes}.
#' @param sampling_parameters List. Passed to the sample_model function.
#' @param model_settings List. Passed to the sample_model function.
#' @param sample_model Function such as \link{sample_log_normal_model}.
#' @param pedigree (optionally) [ped][pedtools::ped] object. Contributors can be named pedigree members.
#' @param results_directory (optionally) Character with path to directory where results are written to disk.
#' @param seed (optionally) Integer seed value that can be used to get reproducible runs. If results are written to disk, the 'Run details.txt' file will contain a seed that can be used for reproducing the result.
#' @param number_of_replicates Integer. Number of replicate simulations for each sample.
#' @param write_non_contributors Logical. If TRUE, sampled genotypes for non-contributing pedigree members will also be written to disk. Defaults to FALSE.
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
#' @examples
#' freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc.csv",
#'                            package = "simDNAmixtures"))
#' gf <- gf_configuration()
#'
#' sampling_parameters <- list(min_mu = 50., max_mu = 5e3,
#'                            min_cv = 0.05, max_cv = 0.35,
#'                            degradation_shape1 = 10, degradation_shape2 = 1)
#'
#' mixtures <- sample_mixtures(n = 2, contributors = c("U1", "U2"), freqs = freqs,
#'                             sampling_parameters = sampling_parameters,
#'                             model_settings = gf$gamma_settings_no_stutter,
#'                             sample_model = sample_gamma_model)
#'
#'# sample a mixture of two siblings taking into account
# linkage between vWA and D12
#'
#' linkage_map <- data.frame(chromosome = c("12","12"),
#'                           locus = c("vWA", "D12391"),
#'                           position = c(16.56662766, 29.48590551))
#'
#' ped_sibs <- pedtools::nuclearPed(children = c("Sib1", "Sib2"))
#'
#' sibs_mix <- sample_mixtures(n = 1, contributors = c("Sib1", "Sib2"),
#'                             freqs = freqs,
#'                             linkage_map = linkage_map,
#'                             pedigree = ped_sibs,
#'                             sampling_parameters = sampling_parameters,
#'                             model_settings = gf$gamma_settings_no_stutter,
#'                             sample_model = sample_gamma_model)
#' @export
sample_mixtures <- function(n, contributors, freqs,
                            linkage_map,
                            sampling_parameters, model_settings,
                            sample_model, pedigree,
                            results_directory,
                            seed,
                            number_of_replicates = 1L,
                            write_non_contributors = FALSE,
                            tag = "simulation", silent = FALSE){

  # validate inputs
  .validate_integer(n, "n", require_strictly_positive = TRUE)
  .validate_logical(write_non_contributors, "write_non_contributors")
  .validate_integer(number_of_replicates, "number_of_replicates",
                    require_strictly_positive = TRUE)

  # if seed is missing we obtain one here
  seed_validated <- .validate_or_generate_seed(seed)
  set.seed(seed_validated)

  number_of_contributors <- length(contributors)

  # init results dir if used
  write_to_disk <- FALSE
  if (!missing(results_directory)){
    write_to_disk <- TRUE

    results_dirs <- .init_results_directory(results_directory, tag, seed_validated, pedigree)
  }

  # pre allocate list for sampling results and sample names
  samples <- vector(mode = "list", length = n * number_of_replicates)
  sample_names <- character(n * number_of_replicates)

  # keep track of sample index even if replicates are used
  i_sample_out <- 1L
  for (i_sample in seq_len(n)){
    all_genotypes <- sample_contributor_genotypes(contributors, freqs,
                                                  linkage_map = linkage_map,
                                                  pedigree = pedigree,
                                                  loci = model_settings$locus_names,
                                                  return_non_contributors = write_non_contributors)

    contributor_genotypes <- all_genotypes[contributors]

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
        if (write_non_contributors){
          write_knowns(all_genotypes, results_dirs$knowns_dir, sample_name)
        }
        else{
          write_knowns(contributor_genotypes, results_dirs$knowns_dir, sample_name)
        }
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
