#' @title Sample mixtures with random genotypes and random parameters according to priors
#'
#' @param n Integer. Number of samples.
#' @param contributors Character vector with unique names of contributors. Valid names are "U1", "U2", ... for unrelated contributors or the names of pedigree members for related contributors.
#' @param freqs Allele frequencies (see \link{read_allele_freqs})
#' @param sampling_parameters List. Passed to the sample_model function.
#' @param model_settings List. Passed to the sample_model function.
#' @param sample_model Function such as \link{sample_log_normal_model}.
#' @param pedigree (optionally) \link[pedtools]{ped} object.
#' @param results_directory (optionally) Character with path to directory where results are written to disk.
#' @param tag Character. Used for sub directory name when results_directory is provided.
#' @examples
#' freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc.csv",package = "SimMixDNA"))
#' gf <- get_GlobalFiler_3500_data()
#'
#' sampling_parameters <- list(min_mu = 50., max_mu = 5e3,
#'                            min_cv = 0.05, max_cv = 0.35,
#'                            degradation_shape1 = 10, degradation_shape2 = 1)
#'
#' mixtures <- sample_mixtures_with_relatives(n = 10, contributors = c("U1", "U2"), freqs = freqs,
#'                                           sampling_parameters = sampling_parameters,
#'                                           model_settings = gf$gamma_settings_no_stutter,
#'                                           sample_model = sample_gamma_model)
#'
#' @export
sample_mixtures_with_relatives <- function(n, contributors, freqs,
                                           sampling_parameters, model_settings,
                                           sample_model, pedigree,
                                           results_directory,
                                           tag = "simulation"){

  if (length(n) != 1){
    stop("n needs to have length 1")
  }

  if (!(is.numeric(n) | is.integer(n))){
    stop("n needs to be integer valued")
  }

  if (as.character(n) != as.character(as.integer(n))){
    stop("n needs to be integer valued")
  }

  number_of_contributors <- length(contributors)

  write_to_disk <- FALSE
  if (!missing(results_directory)){
    if (!dir.exists(results_directory)){
      dir.create(results_directory, recursive = TRUE)
    }

    write_to_disk <- TRUE

    sub_dir <- file.path(results_directory, paste0(format(Sys.time(), "%Y-%m-%d %H_%M"), " ", tag))

    dir.create(sub_dir, recursive = TRUE)

    mixtures_csv_dir <- file.path(sub_dir,"mixtures_csv")
    dir.create(mixtures_csv_dir,recursive = TRUE)

    annotated_mixtures_dir <- file.path(sub_dir,"annotated_mixtures")
    dir.create(annotated_mixtures_dir,recursive = TRUE)

    knowns_dir <- file.path(sub_dir,"knowns")
    dir.create(knowns_dir,recursive = TRUE)
  }

  samples <- list()

  sample_names <- character(n)

  for (i_sample in seq_len(n)){

    sample_name <- paste0("sample", "_", sprintf("%04d", i_sample))

    sample_names[i_sample] <- sample_name

    contributor_genotypes <- sample_contributor_genotypes(contributors, freqs, pedigree,
                                                          loci = model_settings$locus_names)

    model <- sample_model(number_of_contributors = number_of_contributors,
                          sampling_parameters = sampling_parameters,
                          model_settings = model_settings)

    annotated_mixture <- sample_mixture_from_genotypes(contributor_genotypes, model, sample_name)

    mixture <- get_bare_mixture(annotated_mixture)

    samples[[i_sample]] <- list(sample_name = sample_name,
                                contributor_genotypes = contributor_genotypes,
                                model = model,
                                annotated_mixture = annotated_mixture,
                                mixture = mixture)

    if (write_to_disk){
      annnotated_path <- file.path(annotated_mixtures_dir, paste0(sample_name,"_annotated.csv"))
      write.csv(x = annotated_mixture, file = annnotated_path,
                quote = FALSE, row.names = FALSE, na = "")

      mixture_path <- file.path(mixtures_csv_dir, paste0(sample_name,".csv"))
      write.csv(x = mixture, file = mixture_path,
                quote = FALSE, row.names = FALSE, na = "")

      write_knowns(contributor_genotypes, knowns_dir, sample_name)

      db_path <- file.path(sub_dir, "references.csv")
      write_knowns_as_reference_db(samples, db_path)
    }
  }

  names(samples) <- sample_names

  # prepare additional outputs
  parameter_summary <- get_parameter_summary(samples)
  smash <- get_SMASH_from_samples(samples)
  table <- SMASH_to_wide_table(smash)

  if (write_to_disk){
    parameter_summary_path <- file.path(sub_dir, "parameter_summary.csv")
    write.csv(parameter_summary, file = parameter_summary_path,
              quote = FALSE, na = "", row.names = FALSE)

    smash_path <- file.path(sub_dir, paste0(tag, " SMASH.csv"))
    write.csv(smash, file = smash_path,
              quote = FALSE, na = "", row.names = FALSE)

    table_path <- file.path(sub_dir, paste0(tag, " table.txt"))
    write.table(x = table,
                file = table_path, quote = FALSE,
                sep = "\t", row.names = FALSE, na = "")
  }

  list(call = match.call(),
       samples = samples,
       smash = smash,
       parameter_summary = parameter_summary)
}
