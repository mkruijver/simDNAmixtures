# used by the sample_mixtures and sample_mixtures_from_genotypes functions
# to store "user friendly" outputs generated from simulation results
.init_results_directory <- function(results_directory, tag, seed_validated, pedigree){

  if (!dir.exists(results_directory)){
    dir.create(results_directory, recursive = TRUE)
  }

  sub_dir <- file.path(results_directory, paste0(
    format(Sys.time(), "%Y-%m-%d %H_%M_%S"), " ", tag))

  dir.create(sub_dir, recursive = TRUE)

  run_details_file <- file.path(sub_dir, "Run Info.txt");
  write(c(paste0("Simulation started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
          "Call: ",
          deparse(match.call()),"",
          paste0("Seed: ", seed_validated)
  ),
  file = run_details_file)

  if (!missing(pedigree)){
    .write_pedigree(pedigree, file.path(sub_dir,"Pedigree.csv"))
  }

  mixtures_csv_dir <- file.path(sub_dir,"Mixtures csv")
  dir.create(mixtures_csv_dir,recursive = TRUE)

  mixtures_wide_dir <- file.path(sub_dir,"Mixtures table")
  dir.create(mixtures_wide_dir,recursive = TRUE)

  annotated_mixtures_dir <- file.path(sub_dir,"Mixtures annotated")
  dir.create(annotated_mixtures_dir,recursive = TRUE)

  knowns_dir <- file.path(sub_dir,"References by mixture")
  dir.create(knowns_dir,recursive = TRUE)

  list(sub_dir = sub_dir,
       run_details_file = run_details_file,
       mixtures_csv_dir = mixtures_csv_dir,
       mixtures_wide_dir = mixtures_wide_dir,
       annotated_mixtures_dir = annotated_mixtures_dir,
       knowns_dir = knowns_dir)
}

.write_pedigree <- function(pedigree, filename){
  utils::write.csv(as.data.frame(pedigree),
                   file = filename,
                   quote = FALSE, row.names = FALSE)
}

.write_mixture <- function(sample, results_dirs){

  sample_name <- sample$sample_name

  ## annotated
  annotated_path <- file.path(results_dirs$annotated_mixtures_dir, paste0(sample_name," annotated.csv"))
  utils::write.csv(x = sample$annotated_mixture, file = annotated_path,
                   quote = FALSE, row.names = FALSE, na = "")

  ## csv
  mixture_csv_path<- file.path(results_dirs$mixtures_csv_dir, paste0(sample_name,".csv"))
  utils::write.csv(x = sample$mixture, file = mixture_csv_path,
                   quote = FALSE, row.names = FALSE, na = "")

  smash_sample <- get_SMASH_from_samples(list(sample))
  table_sample <- SMASH_to_wide_table(smash_sample)

  ## txt (wide table)
  mixture_wide_path <- file.path(results_dirs$mixtures_wide_dir, paste0(sample_name,".txt"))
  utils::write.table(x = table_sample,
                     file = mixture_wide_path, quote = FALSE,
                     sep = "\t", row.names = FALSE, na = "")
}

.prepare_summaries <- function(samples){
  parameter_summary <- get_parameter_summary(samples)
  smash <- get_SMASH_from_samples(samples)
  table <- SMASH_to_wide_table(smash)
  db <- get_knowns_as_reference_db(samples)

  list(parameter_summary = parameter_summary,
       smash = smash,
       table = table,
       db = db)
}

.write_summaries <- function(summaries, tag, results_dirs){
  parameter_summary_path <- file.path(results_dirs$sub_dir, "Parameter Summary.csv")
  utils::write.csv(summaries$parameter_summary, file = parameter_summary_path,
                   quote = FALSE, na = "", row.names = FALSE)

  smash_path <- file.path(results_dirs$sub_dir, paste0(tag, " SMASH.csv"))
  utils::write.csv(summaries$smash, file = smash_path,
                   quote = FALSE, na = "", row.names = FALSE)

  table_path <- file.path(results_dirs$sub_dir, paste0(tag, " table.txt"))
  utils::write.table(x = summaries$table,
                     file = table_path, quote = FALSE,
                     sep = "\t", row.names = FALSE, na = "")

  ## all knowns as single db (csv)
  db_path <- file.path(results_dirs$sub_dir, "References DB.csv")

  utils::write.csv(x = summaries$db, file = db_path,
                   row.names = FALSE, quote = FALSE)
}
