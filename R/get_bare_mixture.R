get_bare_mixture <- function(annotated_mixture){
  idx <- annotated_mixture$HeightAtOrAboveDetectionThreshold

  data.frame(Locus = annotated_mixture$Marker[idx],
             Allele = annotated_mixture$Allele[idx],
             Height = as.integer(round(annotated_mixture$Height[idx],digits = 0)),
             Size = round(annotated_mixture$Size[idx],2))
}
