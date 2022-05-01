#' Converts SMASH (SampleName, Marker, Allele, Size, Height) data to a wide table
#'
#' @param x DataFrame with SampleName, Marker, Allele, Size, Height columns
#' @return DataFrame with columns: Sample Name, Marker, Allele 1, Allele 2, ..., Size 1, Size 2, ..., Height 1, Height 2, ...
#' @examples
#' # generate example data in SMASH form
#' freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc.csv",
#' package = "simDNAmixtures"))
#' data(gf)
#'
#' sampling_parameters <- list(min_mu = 50., max_mu = 5e3,
#'                             min_cv = 0.05, max_cv = 0.35,
#'                             degradation_shape1 = 10, degradation_shape2 = 1)
#'
#' mixtures <- sample_mixtures(n = 2, contributors = c("U1"), freqs = freqs,
#'                             sampling_parameters = sampling_parameters,
#'                             model_settings = gf$gamma_settings_no_stutter,
#'                             sample_model = sample_gamma_model)
#'
#' # convert from SMASH to wide table
#' wide_table <- SMASH_to_wide_table(mixtures$smash)
#' @export
SMASH_to_wide_table <- function(x) {

  if (!is.data.frame(x)){
    stop("x is not a DataFrame")
  }

  mandatory_columns <- c("SampleName", "Marker", "Allele", "Size", "Height")
  for (col in mandatory_columns){
    if (!(col %in% names(x))){
      stop("x does not have a ", col, " column")
    }
  }

  # determine the maximum number of alleles at a marker across samples
  x_by_sample_name <- split(x, x$SampleName)
  max_number_of_alleles <- max(sapply(x_by_sample_name, function(y) max(table(y$Marker))))

  wide <- get_empty_wide_table(0, max_number_of_alleles)

  for (i_sample in seq_along(x_by_sample_name)){

    x_sample <- x_by_sample_name[[i_sample]]

    sample_name <- names(x_by_sample_name)[i_sample]

    markers <- unique(x_sample$Marker)

    x_sample_wide <- get_empty_wide_table(length(markers), max_number_of_alleles)

    x_sample_wide[["Sample Name"]] <- sample_name
    x_sample_wide$Marker <- markers

    # ensure markers are in order
    x_sample <- x_sample[order(match(x_sample$Marker, markers)),]

    # determine allele number by row number
    allele_number_by_row_idx <- sequence(rle(x_sample$Marker)$lengths)

    for (i_row in seq_len(nrow(x_sample))){
      i_wide_row <- match(x_sample$Marker[i_row], markers)

      allele_number <- allele_number_by_row_idx[i_row]

      x_sample_wide[i_wide_row, 2 + allele_number] <- x_sample$Allele[i_row]
      x_sample_wide[i_wide_row, 2 + allele_number + max_number_of_alleles] <- x_sample$Size[i_row]
      x_sample_wide[i_wide_row, 2 + allele_number + 2 * max_number_of_alleles] <- x_sample$Height[i_row]
    }

    wide <- rbind(wide, x_sample_wide)
  }

  wide
}

get_empty_wide_table <- function(number_of_rows, number_of_alleles){

  empty_character_column <- rep(NA_character_, number_of_rows)
  empty_numeric_column <- rep(NA_real_, number_of_rows)

  wide <- data.frame("Sample Name" = empty_character_column,
                     Marker = empty_character_column,
                     stringsAsFactors = FALSE, check.names = FALSE)

  # Allele1, Allele2, ...
  for (i_allele in seq_len(number_of_alleles)){
    wide[[paste0("Allele", i_allele)]] <- empty_character_column
  }

  # Size1, Size2, ...
  for (i_allele in seq_len(number_of_alleles)){
    wide[[paste0("Size", i_allele)]] <- empty_numeric_column
  }

  # Height1, Height2, ...
  for (i_allele in seq_len(number_of_alleles)){
    wide[[paste0("Height", i_allele)]] <- empty_numeric_column
  }

  wide
}
