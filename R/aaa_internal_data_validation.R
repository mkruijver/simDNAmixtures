#' Validate Allele Frequencies
#'
#' This is an internal function to check the format and contents of the frequency list.
#'
#' @param freqs List of allele frequencies.
#' @param required_loci Character vector of required loci names.
#' @keywords internal
.validate_freqs <- function(freqs, required_loci = NULL){
  if (!is.list(freqs)){
    stop("freqs should be a list")
  }
  if (!all(sapply(freqs, is.numeric))){
    stop("freqs should be a list of numeric vectors")
  }

  if (!is.null(required_loci)){
    for (locus in required_loci){
      if (!(locus %in% names(freqs))){
        stop(paste0("freqs not available for locus "), locus)
      }
    }
  }
}

.validate_sex_locus_freqs <- function(a, f){

  sum_of_freqs <- sum(f)
  if (abs(sum_of_freqs - 1) > 1e-5){
    stop("Sum of freqs at sex locus needs to be 1 but is ", sum_of_freqs)
  }

  if (!all(2 == sapply(a, function(x) length(strsplit(x, split = ",")[[1]])))){
    stop("Sex locus needs to have frequency for genotypes of exactly two alleles")
  }

  if (!(AMEL_XX_packed %in% a)){
    stop("Sex locus needs to have frequency for X,X")
  }
  if (!(AMEL_XY_packed %in% a)){
    stop("Sex locus needs to have frequency for X,Y")
  }

  sum_of_xx_xy_freqs <- as.numeric(f[AMEL_XX_packed] + f[AMEL_XY_packed])
  if (abs(sum_of_xx_xy_freqs - 1) > 1e-5){
    stop("Sum of X,X and X,Y freqs at sex locus needs to be 1 but is ", sum_of_xx_xy_freqs)
  }
}

.ped_forbidden_names <- paste0("U", seq_len(10))

#' Validate Pedigree Object
#'
#' This internal function checks if the input object is of class `ped`. If not, it stops execution with an error message. Additionally, it can check for forbidden names in the pedigree IDs if specified.
#'
#' @param pedigree An object that is expected to be of class `ped`.
#' @param disallow_U_names A logical flag indicating whether to disallow IDs starting with "U" followed by a digit (e.g., "U1", "U2", ..., "U10") in the pedigree. Default is \code{FALSE}.
#' @return This function does not return a value. It stops execution if the input is not of class `ped` or if forbidden names are found when \code{disallow_U_names} is \code{TRUE}.
#' @keywords internal
#' @examples
#' \dontrun{
#'   pedigree <- create_ped() # Assuming create_ped() returns an object of class 'ped'
#'   .validate_pedigree(pedigree)
#'   .validate_pedigree(pedigree, disallow_U_names = TRUE)
#' }
.validate_pedigree <- function(pedigree, disallow_U_names = FALSE){

  if (!inherits(pedigree, "ped")){
    stop("pedigree should be of class ped")
  }

  if (disallow_U_names){
    ped_names <- pedigree$ID

    forbidden_names_in_ped <- intersect(ped_names, .ped_forbidden_names)

    if (length(forbidden_names_in_ped) > 0){
      stop("Pedigree contains illegal name(s): ", paste(forbidden_names_in_ped, collapse = ", "))
    }

  }

}

.validate_linkage_map <- function(linkage_map){#, loci = NULL){

  if (!is.data.frame(linkage_map)){
    stop("linkage_map should be a DataFrame")
  }

  required_columns <- c("chromosome", "locus", "position")
  for (required_column in required_columns){
    x <- linkage_map[[required_column]]
    if (is.null(x)){
      stop("linkage_map should be a DataFrame with a column named: ", required_column)
    }
  }
}

.validate_logical <- function(x, param_name){
  if (missing(param_name)) param_name <- as.character(match.call()$x)

  if (!is.logical(x)){
    stop(param_name, " needs to be a logical")
  }

  if (length(x) != 1){
    stop(param_name, " needs to be a logical of length 1")
  }

  x
}

.validate_integer <- function(n, param_name, require_strictly_positive = FALSE){
  if (missing(param_name)) param_name <- as.character(match.call()$n)

  if (length(n) != 1){
    stop(param_name, " needs to have length 1")
  }

  if (!(is.numeric(n) | is.integer(n))){
    stop(param_name, " needs to be integer valued")
  }

  if (as.character(n) != as.character(as.integer(n))){
    stop(param_name, " needs to be integer valued")
  }

  if (require_strictly_positive & (n <= 0)){
    stop(param_name, " needs to be strictly positive")
  }

  as.integer(n)
}

.validate_character <- function(x, param_name, required_length){
  if (missing(param_name)) param_name <- as.character(match.call()$x)

  if (!missing(required_length)){
    if (length(x) != required_length){
      stop(param_name, " needs to have length ", required_length)
    }
  }

  if (!is.character(x)){
    stop(param_name, " needs to be a character")
  }
}

.validate_numeric <- function(x, param_name, require_strictly_positive = FALSE,
                              require_nonnegative = FALSE,
                              require_length_one = TRUE){

  if (missing(param_name)) param_name <- as.character(match.call()$x)

  if (require_length_one & (length(x) != 1)){
    stop(param_name, " needs to have length 1")
  }

  if (!(is.numeric(x) | is.integer(x))){
    stop(param_name, " needs to be numeric")
  }

  if (require_nonnegative & any(x < 0.)){
    stop(param_name, " needs to be non-negative")
  }

  if (require_strictly_positive & any(x <= 0.)){
    stop(param_name, " needs to be strictly positive")
  }

  as.numeric(x)
}

.validate_clamp <- function(x, param_name, minimum = NA, maximum = NA,
                            strict_minimum = FALSE,
                            strict_maximum = FALSE){

  if (missing(param_name)) param_name <- as.character(match.call()$x)

  if (!is.na(minimum)){
    if (any(x < minimum)){
      stop(param_name, " needs to be at least ", minimum)
    }
    if (strict_minimum & any(x <= minimum)){
      stop(param_name, " needs to be more than ", minimum)
    }
  }

  if (!is.na(maximum)){
    if (any(x > maximum)){
      stop(param_name, " needs to be at most ", maximum)
    }
    if (strict_maximum & any(x >= strict_maximum)){
      stop(param_name, " needs to be less than ", strict_maximum)
    }
  }
}

.validate_or_generate_seed <- function(seed){
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
    return(seed)
  }
  else{

    # pick a seed up to 1 million
    # and return this seed for reproducible results even if no seed provided
    seed <- sample.int(n = 1e6, size = 1)

    return(seed)
  }
}

.parse_STRmix_boolean <- function(x){
  isTRUE(x[[1]] == "Y")
}

.parse_STRmix_double <- function(x){
  as.numeric(strsplit(x, split = ",")[[1]])
}


#
#   if (!is.null(loci)){
#     missing_loci <- loci[!loci %in% linkage_map$locus]
#
#     stop("loci missing from linkage map:", missing_loci)
#   }
