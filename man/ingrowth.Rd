% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ingrowth.R
\name{ingrowth}
\alias{ingrowth}
\title{Ingrowth imputation for permanent sample plots}
\usage{
ingrowth(treedata, plotdata, plot, meas, tree, diam, thresh)
}
\arguments{
\item{treedata}{Tree data frame. Contains at least plot, measurement,
and tree identifiers, and diameter values.}

\item{plotdata}{Plot data frame. Contains at least a plot identifier,
and the minimum diameter threshold for the plot.}

\item{plot}{Name of the plot identifier in \code{treedata} and
\code{plotdata}.}

\item{meas}{Name of the measurement identifier in \code{treedata}.}

\item{tree}{Name of the tree identifier in \code{treedata}.}

\item{diam}{Name of the diameter values in \code{treedata}.}

\item{thresh}{Name of the minimum diameter threshold in \code{treedata}.}
}
\value{
A data frame with \code{treedata} followed by the imputed ingrowth
tree records. A logical variable named \code{ingrowth} is appended,
containing \code{FALSE} for the observed trees and \code{TRUE}
for the ingrowth. The values in the ingrowth records for
variables other than the identifiers and diameters are inherited
from later measurements and are often meaningless.
}
\description{
Adds ingrowth records to a tree data frame
}
\details{
The data should not include plots with substantial ingrowth at
the oldest measurement. The measurement identifier values must
be capable of being sorted in time, e.g., they could be
consecutive integer indices, years, or dates. The other
identifiers can be of any type (numbers, strings, or factors).
It is advisable to discard measurements where the imputed
ingrowth basal area is a large fraction of the total.
}
\examples{
 withIngrowth <- ingrowth(treeData, plotData, plot="PlotId",
                          meas = "MeasId", tree = "TreeId",
                          diam = "DBH", thresh = "Thresh")

}
\references{
García, O. (2021)"Imputing Ingrowth in Even-aged Permanent
Sample Plots". \emph{Mathematical and Computational Forestry
& Natural-Resource Sciences (MCFNS)} (to appear)
}
