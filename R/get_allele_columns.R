
.get_allele_columns <- function(y){

  y[startsWith(names(y), "Allele")]
}
