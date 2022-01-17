
get_allele_index <- function(x, marker, allele){
  which(x$Marker == marker & x$Allele == allele)
}

get_stutter_target <- function(parent, delta){
  as.character(as.numeric(parent) + delta)
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

set_or_add_df_variable <- function(x, marker, allele, size, value, column_name){
  idx <- get_allele_index(x, marker, allele)

  if(length(idx)==0){

    new_df <- data.frame(Marker=marker, Allele=allele,
                         Size=size, stringsAsFactors = FALSE)
    new_df[[column_name]] <- value

    return(dplyr::bind_rows(x, new_df))
  }else if(length(idx)==1){
    x[[column_name]][idx] <- value
    return(x)
  }else{
    stop("something wrong")
  }
}

