% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/precipAggregate.R
\name{precip_hour2daily}
\alias{precip_hour2daily}
\title{Precipitation aggregation from hourly to daily}
\usage{
precip_hour2daily(
  start_hour,
  end_hour,
  dirHour,
  inFormat,
  dirOUT,
  outFormat,
  obs_hour = 8,
  min_frac = 0.9,
  time_zone = "Africa/Kigali",
  long_name = "Merged AWS-Radar daily rainfall"
)
}
\arguments{
\item{start_hour, end_hour}{start and end time of the period to aggregate 
same time zone as \code{time_zone}, format "YYYY-mm-dd HH:00"}

\item{dirHour}{full path to the directory containing the NetCDF files}

\item{inFormat}{format of the input file name. Ex: "precip_adj_\%s\%s\%s\%s.nc"}

\item{dirOUT}{full path to the directory to save the aggregate data}

\item{outFormat}{format of the output file name. Ex: "precip_adj_\%s\%s\%s.nc"}

\item{obs_hour}{observation hour in local time. Default 8h for Rwanda}

\item{min_frac}{minimum fraction of non missing values}

\item{time_zone}{the time zone of \code{start_time}, \code{end_time} and the input hourly netCDF files.
Options: "Africa/Kigali" or "UTC". Default "Africa/Kigali"}

\item{long_name}{long name for the NetCDF variable precip}
}
\description{
Aggregate hourly precipitation to daily.
}
