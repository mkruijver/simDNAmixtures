# simDNAmixtures 1.1.1 (April 15, 2025)
* Require R >= 4.1.0 to address NOTE

# simDNAmixtures 1.1.0 (April 14, 2025)
* Add `sample_many_pedigree_genotypes` function to more efficiently sample many genotypes
* Mixtures can now be simulated including the AMEL locus
* YSTR loci can now be simulated (see vignette)
* `sample_mixtures` now supports replicate simulations through `number_of_replicates` parameter
* `sample_contributors` now supports linkage
* The drop model model (`drop_model`) can now be used instead of `gamma_model` or `log_normal_model`
* Mixture proportions/template parameters are no longer sorted in decreasing order by `sample_gamma_model` and `sample_log_normal_model`
* Replace data(gf) to runtime gf_configuration() function to fix deserialised functions throwing warnings
* Add exemplar data for Kintelligence (map and 1000G allele frequencies)
* Bug fixed in `read_allele_freqs` where `N` attribute was removed if `normalise = TRUE` or `remove_zeroes = TRUE`

# simDNAmixtures 1.0.1 (May 3, 2022)
* Removed `\dontrun` from `gamma_model` example

# simDNAmixtures 1.0 (May 1, 2022)

* Initial CRAN submission
* Added a `NEWS.md` file to track changes to the package
