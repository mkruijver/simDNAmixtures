get_SMASH_from_samples <- function(samples){

  smash_by_sample <- lapply(samples, function(s){
    data.frame(SampleName = rep(s$sample_name, nrow(s$mixture)),
               Marker = s$mixture$Locus,
               Allele = s$mixture$Allele,
               Size = s$mixture$Size,
               Height = s$mixture$Height)
  })

  smash <- do.call(rbind, smash_by_sample)
  rownames(smash) <- NULL

  smash
}
