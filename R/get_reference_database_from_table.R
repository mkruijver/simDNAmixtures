# makes a reference database from collated knowns
get_reference_database <- function(x){

  x_by_profile <- split(x, f =  paste0(x$CaseNumber,"!",x$`Sample Name`))

  y <- x_by_profile[[1]]

  db_rows <- lapply(x_by_profile, function (y){
    data.frame(CaseNumber = y$CaseNumber[1], "Sample Name" = y$`Sample Name`[1] ,
               data.frame(t(stats::setNames(as.vector(rbind(y$Allele1, y$Allele2)),
                                     paste0(rep(y$Locus, each=2),c("","____2____")))
                            # rep(y$Locus, each=2)
               ),
               check.names = FALSE, stringsAsFactors = FALSE),
               stringsAsFactors = FALSE, check.names = FALSE)
  })

  db <- dplyr::bind_rows(db_rows)
  names(db) <- gsub("____2____", replacement = "", names(db))

  db
}
