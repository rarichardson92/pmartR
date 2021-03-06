% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trelliVis.R
\name{set_ylimits}
\alias{set_ylimits}
\title{Sets y-axis limits for trellData plotting}
\usage{
set_ylimits(yvalues, increment, y_range = NULL, y_max = NULL,
  y_min = NULL, include_zero = TRUE)
}
\arguments{
\item{yvalues}{y-values for plotting}

\item{increment}{An increment set based on the maximum/minimum y-values}

\item{y_range}{Specify a range for the plot y-axis. Will calculate the range based on one of y_max or y_min parameters or from the median of y-values where y_max and y_min are not defined.}

\item{y_max}{Sets the maximum y-value for the y-axis.}

\item{y_min}{Sets the minimum y-value for the y-axis.}

\item{include_zero}{Should zero be included regardless of scaling?}
}
\description{
Sets y-axis limits for trellData plotting. Used by plot_comp.
}
\seealso{
\link[pmartR]{format_plot}

\link[pmartR]{list_y_limits}

\link[pmartR]{set_increment}
}
\author{
Rachel Richardson
}
