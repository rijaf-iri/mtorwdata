% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validationStats.R
\name{stats_volumetric_quantile}
\alias{stats_volumetric_quantile}
\title{Compute volumetric statistics.}
\usage{
stats_volumetric_quantile(x, y, thres)
}
\arguments{
\item{x}{matrix of the reference observations data used to validate the simulated or estimated data.\cr
The row represents the observations \cr
The column represents the points or stations}

\item{y}{matrix, same dimension as \code{x}, containing the simulated or estimated data to validate.}

\item{thres}{vector of thresholds to use for each stations. \cr
If one value is provided, all the stations will use this threshold \cr
If a vector is provided, it must be the same length of the column of \code{x}}
}
\value{
A list object
\itemize{
\item{\strong{statistics}: }{a matrix containing the values of the computed statistics.\cr
                    The column represents the points/stations and the row for the statistics}
\item{\strong{description}: }{vector describing the statistics}
\item{\strong{perfect.score}: }{vector of the perfect score values for each statistics}
}
}
\description{
Compute volumetric statistics using a quantile thresholds for validation.
}
