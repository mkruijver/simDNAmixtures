get_allele_index <- function(x, marker, allele){
  which(x$Marker == marker & x$Allele == allele)
}

repeats_to_bp <- function(repeats, repeat_length){
  whole_repeats <- floor(repeats)
  partials <- 10.0 * (repeats - whole_repeats)

  whole_repeats * repeat_length + partials
}

bp_to_repeats <- function(bp, repeat_length){
  if (abs(round(bp)-bp) > 1e-9) stop("Expected integer number of basepairs")
  bp <- round(bp)

  partials <- bp %% repeat_length
  whole_repeats <- (bp - partials) / repeat_length

  whole_repeats + 0.1 * partials
}

# converts e.g. 9.3 to 9.74 for a tetranucleotide
repeats_to_decimals <- function(repeats, repeat_length){
  bp <- repeats_to_bp(repeats, repeat_length)

  bp / repeat_length
}

get_stutter_target <- function(parent, delta, repeat_length){

  if (length(delta) == 1){
    return(as.character(as.numeric(parent) + delta))
  }
  else if(length(delta) == 2){
    bp_parent <- repeats_to_bp(as.numeric(parent), repeat_length)
    bp_target <- bp_parent + delta[1] * repeat_length + delta[2]

    target <- as.character(bp_to_repeats(bp_target, repeat_length))

    return(target)
  }
  else{
    stop("delta is not length one or two")
  }
}

add_expected_allelic_peak_height <- function(x, marker, allele, size, expected){
  idx <- get_allele_index(x, marker, allele)

  if(length(idx)==0){
    return(dplyr::bind_rows(x, data.frame(Marker=marker, Allele=allele,
                                          Size=size,
                                          ExpectedAllele=expected, stringsAsFactors = FALSE)))
  }else if(length(idx)==1){
    x$ExpectedAllele[idx] <- x$ExpectedAllele[idx] + expected
    return(x)
  }else{
    stop("something wrong")
  }
}

add_expected_peak_height <- function(x, marker, allele, size, expected, column_name){
  idx <- get_allele_index(x, marker, allele)

  if(length(idx)==0){

    new_df <- data.frame(Marker=marker, Allele=allele,
                         ExpectedAllelePreStutter = 0.,
                         ExpectedAllele = 0., ExpectedStutter = 0.,
                         Size=size, stringsAsFactors = FALSE)
    new_df[[column_name]] <- expected

    return(dplyr::bind_rows(x, new_df))
  }else if(length(idx)==1){
    x[[column_name]][idx] <- sum(x[[column_name]][idx], expected, na.rm = TRUE)
    return(x)
  }else{
    stop("something wrong")
  }
}

set_or_add_df_variable <- function(x, marker, allele, size, value, column_name,
                                   sum = FALSE){
  idx <- get_allele_index(x, marker, allele)

  if(length(idx)==0){

    new_df <- data.frame(Marker=marker, Allele=allele,
                         Size=size, stringsAsFactors = FALSE)
    new_df[[column_name]] <- value

    return(dplyr::bind_rows(x, new_df))
  }else if(length(idx)==1){

    if (!sum){
      x[[column_name]][idx] <- value
    }
    else{
      x[[column_name]][idx] <- x[[column_name]][idx] + value
    }
    return(x)
  }else{
    stop("something wrong")
  }
}

.extract_repeat_length_by_locus_from_STRmix_kit <- function(kit_xml){
  repeat_lengths <- sapply(kit_xml$profilingKit$loci, function(x) as.numeric(x$repeatLength[[1]]))

  # in a different version of the kit xml the repeat length is found as an attribute
  if (length(unlist(repeat_lengths)) == 0){
    repeat_lengths <- as.numeric(unlist(sapply(kit_xml$profilingKit$loci, function(x) {
      repeat_length <- attr(x, "repeatLength")
      if (is.null(repeat_length)) 0. else repeat_length
    })))
  }
  repeat_length_by_locus <- stats::setNames(repeat_lengths, locus_names)
}
