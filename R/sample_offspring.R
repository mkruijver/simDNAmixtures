#' @title Sample offspring from two parental genotypes
#'
#' @param father DataFrame (see \link{sample_genotype})
#' @param mother DataFrame (see \link{sample_genotype})
#' @param label Character SampleName of child
#' @details A genotype is sampled according to Mendelian inheritance.
#' @examples
#' # below we read an allele freqs and sample a genotype
#' filename <- system.file("extdata","FBI_extended_Cauc.csv",
#'                         package = "simDNAmixtures")
#' freqs <- read_allele_freqs(filename)
#'
#' # sample parents
#' father <- sample_genotype(freqs, loci = c("D3S1358", "vWA"))
#' mother <- sample_genotype(freqs, loci = c("D3S1358", "vWA"))
#'
#' # sample child
#' child <- sample_offspring(father, mother)
#' @export
sample_offspring <- function(father, mother, label = "Child"){

  if (!is.character(label)){
    stop("label should be a character vector")
  }
  if (length(label)!=1){
    stop("label should be a character vector of length 1")
  }

  check_genotype_df(father, "father")
  check_genotype_df(mother, "mother")

  if (!identical(father$Locus, mother$Locus)){
    stop("father and mother should have identical locus columns")
  }

  loci <- father$Locus
  number_of_loci <- length(loci)

  a <- ifelse(sample(c(FALSE,TRUE), replace = TRUE, size = nrow(father)),father$Allele1, father$Allele2)
  b <- ifelse(sample(c(FALSE,TRUE), replace = TRUE, size = nrow(mother)),mother$Allele1, mother$Allele2)

  # sort alleles
  for (i_locus in seq_len(number_of_loci)){
    a_i <- a[[i_locus]]
    b_i <- b[[i_locus]]

    a_i_numeric <- as.numeric(a_i)
    b_i_numeric <- as.numeric(b_i)

    a_i_is_numeric <- as.character(a_i_numeric) == a_i
    b_i_is_numeric <- as.character(b_i_numeric) == b_i

    if (a_i_is_numeric && b_i_is_numeric && (a_i_numeric > b_i_numeric)){
      a[i_locus] <- b_i
      b[i_locus] <- a_i
    }

  }

  epg <- data.frame("Sample Name" = rep(label, number_of_loci),
                    Locus = loci,
                    Allele1 = a,
                    Allele2 = b,
                    stringsAsFactors = FALSE, check.names = FALSE)

  epg
}

check_genotype_df <- function(x, name){
  if (!is.data.frame(x)){
    stop(paste0(name, " should be DataFrame"))
  }

  # check that all columns are present
  mandatory_columns <- c("Sample Name", "Locus", "Allele1", "Allele2")

  for (column_name in mandatory_columns){
    if (!(column_name %in% names(x))){
      stop(paste0(name, " should have a ", column_name, "column"))
    }
  }

  # check that there is a single sample name
  if (nrow(x) == 0){
    stop(paste0(name, " should have a at least one row"))
  }

  if (!all(x[["Sample Name"]]  == x[["Sample Name"]][1])){
    stop(paste0(name, " should have a unique sample name"))
  }
}
