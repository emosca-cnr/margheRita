% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_NA.R
\name{filter_NA}
\alias{filter_NA}
\title{filter NA}
\usage{
filter_NA(
  mRList = NULL,
  min_metab_in_sample = 100,
  min_sample_with_metab = 3,
  na_value = "NA"
)
}
\arguments{
\item{mRList}{mRList object}

\item{min_metab_in_sample}{min number of metabolites in a sample}

\item{min_sample_with_metab}{min number of samples in which a metabolite must appear}

\item{na_value}{value that indicate missing values}
}
\value{
filtered mRList object
}
\description{
filter NA
}
