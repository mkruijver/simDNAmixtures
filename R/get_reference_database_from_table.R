# makes a reference database from collated knowns
get_reference_database <- function(x){

  x_by_profile <- split(x, f =  paste0(x$CaseNumber,"!",x$`Sample Name`))

  .allele_tables_to_wide_references(x_by_profile)
}
