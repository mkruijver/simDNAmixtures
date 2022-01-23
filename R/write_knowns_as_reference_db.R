
get_knowns_as_reference_db <- function(samples){

  knowns_collated_by_sample <- lapply(samples, function(s){
    g <- do.call(rbind, s$contributor_genotypes)
    g$CaseNumber <- s$sample_name

    g
  })

  knowns_collated <- do.call(rbind, knowns_collated_by_sample)

  db <- get_reference_database(knowns_collated)

  db
}

write_knowns_as_reference_db <- function(samples, path){

  db <- get_knowns_as_reference_db(samples)

  write.csv(x = db, file = path, row.names = FALSE, quote = FALSE)
}
