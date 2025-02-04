% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RLA.R
\name{RLA}
\alias{RLA}
\title{Relative Log Abudance}
\usage{
RLA(
  mRList = NULL,
  include_QC = FALSE,
  logged = FALSE,
  robust = TRUE,
  pal = NULL,
  col_by = "class",
  ...
)
}
\arguments{
\item{mRList}{mRList object}

\item{include_QC}{whether to include or not the QC samples}

\item{logged}{are the input data on log-scale or not?}

\item{robust}{whether to use the median or not}

\item{pal}{color palette}

\item{col_by}{column of sample_ann to color}

\item{...}{further arguments to boxplot function}
}
\value{
mRList object with element "RLA"
}
\description{
Relative Log Abudance
}
