% Generated manually because roxygen2 is unavailable in the execution environment
\name{visualize_airport_delays}
\alias{visualize_airport_delays}
\title{Visualise Average Flight Delays by Destination Airport}
\usage{
visualize_airport_delays()
}
\value{A \code{ggplot2} object visualising mean arrival delays by airport location.}
\description{
Aggregates the \code{nycflights13::flights} data set with \code{dplyr} and
returns a \code{ggplot2} scatter plot showing the mean arrival delay for each
destination airport as a function of its longitude and latitude.
}
\details{
Point colour encodes the average arrival delay while point size reflects the
number of flights contributing to the estimate. The function requires both the
\pkg{dplyr} and \pkg{nycflights13} packages to be installed and available.
}
\examples{
if (requireNamespace("nycflights13", quietly = TRUE)) {
  plt <- visualize_airport_delays()
  print(plt)
}
}
