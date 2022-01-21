write_knowns <- function(contributor_genotypes, knowns_dir, sample_name){
  contributor_names <- names(contributor_genotypes)

  for (contributor_name in contributor_names){
    knowns_path <- file.path(knowns_dir, paste0(sample_name, "_",
                                                contributor_name, ".txt"))

    write.table(x = contributor_genotypes[[contributor_name]],
                file = knowns_path, quote = FALSE,
                sep = "\t",row.names = FALSE)
  }
}
