% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dls_analysis.R
\name{run_all}
\alias{run_all}
\title{Process DLS and titer training data, train model and predict titer loss for new DLS data}
\usage{
run_all(data, lytic, training_pairs, metric, testing_pairs, intercept = FALSE)
}
\arguments{
\item{data}{A dataframe with DLS data. One row per sample}

\item{lytic}{A dataframe with the titer associated with each treatment.}

\item{training_pairs}{A dataframe with two columns, named 'treatment' and 'control'.}

\item{metric}{A string specifying the DLS metric used. Must be a substring of column names in 'data', excluding the 'sample' column.}

\item{testing_pairs}{Optional. A dataframe with two columns, named 'treatment' and 'control'.}
}
\description{
This is a wrapper function that runs all of the analysis, replicating results from the Phage-ELF shiny app.
}
\examples{
run_all(dls_data, lytic_data, training_pairs, 'intens', testing_pairs)
}
