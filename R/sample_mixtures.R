#' @title Sample mixtures with random genotypes and random parameters according to priors
#'
#' @param n Integer. Number of samples.
#' @param contributors Character vector with unique names of contributors. Valid names are "U1", "U2", ... for unrelated contributors or the names of pedigree members for related contributors.
#' @param freqs Allele frequencies (see \link{read_allele_freqs})
#' @param sampling_parameters List. Passed to the sample_model function.
#' @param model_settings List. Passed to the sample_model function.
#' @param sample_model Function such as \link{sample_log_normal_model}.
#' @param pedigree (optionally) \link[pedtools]{ped} object. Contributors can be named pedigree members.
#' @param results_directory (optionally) Character with path to directory where results are written to disk.
#' @param seed (optionally) Integer seed value that can be used to get reproducible runs. If results are written to disk, the 'Run details.txt' file will contain a seed that can be used for reproducing the result.
#' @param write_non_contributors Logical. If TRUE, sampled genotypes for non-contributing pedigree members will also be written to disk. Defaults to FALSE.
#' @param tag Character. Used for sub directory name when results_directory is provided.
#' @examples
#' freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc.csv",package = "SimMixDNA"))
#' gf <- get_GlobalFiler_3500_data()
#'
#' sampling_parameters <- list(min_mu = 50., max_mu = 5e3,
#'                            min_cv = 0.05, max_cv = 0.35,
#'                            degradation_shape1 = 10, degradation_shape2 = 1)
#'
#' mixtures <- sample_mixtures(n = 10, contributors = c("U1", "U2"), freqs = freqs,
#'                             sampling_parameters = sampling_parameters,
#'                             model_settings = gf$gamma_settings_no_stutter,
#'                             sample_model = sample_gamma_model)
#'
#' @export
sample_mixtures <- function(n, contributors, freqs,
                            sampling_parameters, model_settings,
                            sample_model, pedigree,
                            results_directory,
                            seed,
                            write_non_contributors = FALSE,
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

  if (!is.logical(write_non_contributors)){
    stop("write_non_contributors needs to be a logical")
  }

  if (length(write_non_contributors) != 1){
    stop("write_non_contributors needs to be a logical of length 1")
  }

  if (!missing(seed)){

    if (length(seed) != 1){
      stop("seed needs to have length 1")
    }

    if (!(is.numeric(seed) | is.integer(seed))){
      stop("seed needs to be integer valued")
    }

    if (as.character(seed) != as.character(as.integer(seed))){
      stop("seed needs to be integer valued")
    }

    seed <- as.integer(seed)

    set.seed(seed)
  }
  else{

    # pick a seed up to 1 million
    # and return this seed for reproducible results even if no seed provided
    seed <- sample.int(n = 1e6, size = 1)

    set.seed(seed)
  }

  number_of_contributors <- length(contributors)

  write_to_disk <- FALSE
  if (!missing(results_directory)){
    if (!dir.exists(results_directory)){
      dir.create(results_directory, recursive = TRUE)
    }

    write_to_disk <- TRUE

    sub_dir <- file.path(results_directory, paste0(
      format(Sys.time(), "%Y-%m-%d %H_%M_%S"), " ", tag))

    dir.create(sub_dir, recursive = TRUE)

    run_details_file <- file.path(sub_dir, "Run Info.txt");
    write(c(paste0("Simulation started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                 "Call: ",
                 deparse(match.call()),"",
                 paste0("Seed: ", seed)
                 ),
                 file = run_details_file)


    mixtures_csv_dir <- file.path(sub_dir,"Mixtures csv")
    dir.create(mixtures_csv_dir,recursive = TRUE)

    mixtures_wide_dir <- file.path(sub_dir,"Mixtures table")
    dir.create(mixtures_wide_dir,recursive = TRUE)

    annotated_mixtures_dir <- file.path(sub_dir,"Mixtures annotated")
    dir.create(annotated_mixtures_dir,recursive = TRUE)

    knowns_dir <- file.path(sub_dir,"References by mixture")
    dir.create(knowns_dir,recursive = TRUE)
  }

  samples <- list()

  sample_names <- character(n)

  for (i_sample in seq_len(n)){

    all_genotypes <- sample_contributor_genotypes(contributors, freqs, pedigree,
                                                          loci = model_settings$locus_names,
                                                          return_non_contributors = write_non_contributors)

    contributor_genotypes <- all_genotypes[contributors]

    model <- sample_model(number_of_contributors = number_of_contributors,
                          sampling_parameters = sampling_parameters,
                          model_settings = model_settings)

    sample_name <- paste0("sample", "_", sprintf("%04d", i_sample),
                          "_", model$sample_name_suffix)

    sample_names[i_sample] <- sample_name

    annotated_mixture <- sample_mixture_from_genotypes(contributor_genotypes, model, sample_name)

    mixture <- get_bare_mixture(annotated_mixture)

    samples[[i_sample]] <- list(sample_name = sample_name,
                                contributor_genotypes = contributor_genotypes,
                                model = model,
                                annotated_mixture = annotated_mixture,
                                mixture = mixture)

    if (write_to_disk){
      ## annotated
      annnotated_path <- file.path(annotated_mixtures_dir, paste0(sample_name," annotated.csv"))
      utils::write.csv(x = annotated_mixture, file = annnotated_path,
                quote = FALSE, row.names = FALSE, na = "")

      ## csv
      mixture_csv_path<- file.path(mixtures_csv_dir, paste0(sample_name,".csv"))
      utils::write.csv(x = mixture, file = mixture_csv_path,
                quote = FALSE, row.names = FALSE, na = "")

      smash_sample <- get_SMASH_from_samples(samples[i_sample])
      table_sample <- SMASH_to_wide_table(smash_sample)

      ## txt (wide table)
      mixture_wide_path <- file.path(mixtures_wide_dir, paste0(sample_name,".txt"))
      utils::write.table(x = table_sample,
                  file = mixture_wide_path, quote = FALSE,
                  sep = "\t", row.names = FALSE, na = "")

      ## knowns (wide table)
      if (write_non_contributors){
        write_knowns(all_genotypes, knowns_dir, sample_name)
      }
      else{
        write_knowns(contributor_genotypes, knowns_dir, sample_name)
      }


    }
  }

  names(samples) <- sample_names

  # prepare additional outputs
  parameter_summary <- get_parameter_summary(samples)
  smash <- get_SMASH_from_samples(samples)
  table <- SMASH_to_wide_table(smash)

  if (write_to_disk){
    parameter_summary_path <- file.path(sub_dir, "Parameter Summary.csv")
    utils::write.csv(parameter_summary, file = parameter_summary_path,
              quote = FALSE, na = "", row.names = FALSE)

    smash_path <- file.path(sub_dir, paste0(tag, " SMASH.csv"))
    utils::write.csv(smash, file = smash_path,
              quote = FALSE, na = "", row.names = FALSE)

    table_path <- file.path(sub_dir, paste0(tag, " table.txt"))
    utils::write.table(x = table,
                file = table_path, quote = FALSE,
                sep = "\t", row.names = FALSE, na = "")

    ## all knowns as single db (csv)
    db_path <- file.path(sub_dir, "References DB.csv")
    write_knowns_as_reference_db(samples, db_path)

    write(c(paste0("Simulation finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
          file = run_details_file, append = TRUE)

    cat("Finished sampling. Output written to", sub_dir, "\n")
  }

  list(call = match.call(),
       samples = samples,
       smash = smash,
       parameter_summary = parameter_summary)
}
