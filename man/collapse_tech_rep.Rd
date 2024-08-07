% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collapse_tech_rep.R
\name{collapse_tech_rep}
\alias{collapse_tech_rep}
\title{Collapse technical replicates}
\usage{
collapse_tech_rep(mRList = NULL, remove.QC = TRUE)
}
\arguments{
\item{mRList}{mRList object}

\item{remove.QC}{whether to remove QC samples}
}
\value{
mRList object with only biological replicates
}
\description{
Collapse technical replicates of each biological replicate calculating the mean.
The resulting dataframe and metadata contain only biological replicates.
}
