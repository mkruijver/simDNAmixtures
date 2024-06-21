# makes a reference database from collated knowns
# see write_knowns_as_reference_db.R
get_reference_database <- function(x){

  x_by_profile <- split(x, f =  paste0(x$CaseNumber,"!",x$`Sample Name`))

  case_number <- as.vector(sapply(x_by_profile, function(x) x$CaseNumber[1]))
  sample_name <- as.vector(sapply(x_by_profile, function(x) x$`Sample Name`[1]))

  db <- .allele_tables_to_wide_references(x_by_profile)
  data.frame(CaseNumber = case_number, "Sample Name" = sample_name,
             db, check.names = FALSE)
}
