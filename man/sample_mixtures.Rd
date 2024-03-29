% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_mixtures.R
\name{sample_mixtures}
\alias{sample_mixtures}
\title{Sample mixtures with random genotypes and random parameters according to priors}
\usage{
sample_mixtures(
  n,
  contributors,
  freqs,
  sampling_parameters,
  model_settings,
  sample_model,
  pedigree,
  results_directory,
  seed,
  write_non_contributors = FALSE,
  tag = "simulation"
)
}
\arguments{
\item{n}{Integer. Number of samples.}

\item{contributors}{Character vector with unique names of contributors. Valid names are "U1", "U2", ... for unrelated contributors or the names of pedigree members for related contributors.}

\item{freqs}{Allele frequencies (see \link{read_allele_freqs})}

\item{sampling_parameters}{List. Passed to the sample_model function.}

\item{model_settings}{List. Passed to the sample_model function.}

\item{sample_model}{Function such as \link{sample_log_normal_model}.}

\item{pedigree}{(optionally) \link[pedtools:ped]{ped} object. Contributors can be named pedigree members.}

\item{results_directory}{(optionally) Character with path to directory where results are written to disk.}

\item{seed}{(optionally) Integer seed value that can be used to get reproducible runs. If results are written to disk, the 'Run details.txt' file will contain a seed that can be used for reproducing the result.}

\item{write_non_contributors}{Logical. If TRUE, sampled genotypes for non-contributing pedigree members will also be written to disk. Defaults to FALSE.}

\item{tag}{Character. Used for sub directory name when results_directory is provided.}
}
\value{
If \code{results_directory} is provided, this function has the side effect of writing results to disk.

Return value is a list with simulation results:\itemize{
\item \code{call} matched call
\item \code{smash} DataFrame with all samples in SMASH format (see \link{SMASH_to_wide_table})
\item \code{samples} Detailed results for each sample
\item \code{parameter_summary} DataFrame with parameters for each sample
}
}
\description{
Sample mixtures with random genotypes and random parameters according to priors
}
\examples{
freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc.csv",
                           package = "simDNAmixtures"))
data(gf)

sampling_parameters <- list(min_mu = 50., max_mu = 5e3,
                           min_cv = 0.05, max_cv = 0.35,
                           degradation_shape1 = 10, degradation_shape2 = 1)

mixtures <- sample_mixtures(n = 2, contributors = c("U1", "U2"), freqs = freqs,
                            sampling_parameters = sampling_parameters,
                            model_settings = gf$gamma_settings_no_stutter,
                            sample_model = sample_gamma_model)

}
