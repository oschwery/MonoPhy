\name{MonoPhy-package}
\alias{MonoPhy-package}
\alias{MonoPhy}
\docType{package}
\title{
Allows to explore monophyly (or lack of it) of taxonomic groups in a phylogeny.
}
\description{
Requires rooted resolved phylogeny as input and creates a table of genera, their monophyly-status, which taxa cause problems in monophyly etc. Different information can be extracted from the output and a plot function allows visualization of the results in a number of ways.
}
\details{
\tabular{ll}{
Package: \tab MonoPhy\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-22\cr
License: \tab GPL-3\cr
}
The package allows to explore monophyly in phylogenetic trees in a simple and quick way. AssessMonophyly determines whether taxa in a phylogeny are monophyletic. GetSummaryMonophyly extracts a summary table from the results, GetResultMonophyly extracts the main result table, GetIntruderTaxa and GetIntruderTips extract names of taxa which interfer with the monophyly of another taxon and GetAncNodes a table of MRCA nodes. PlotMonophyly finally allows visualizing the results in a couple of different ways.
}
\author{
Orlando Schwery

Maintainer: Orlando Schwery <oschwery@vols.utk.edu>
}
\seealso{
\code{\link{AssessMonophyly}}, \code{\link{GetAncNodes}}, \code{\link{GetIntruderTaxa}}, \code{\link{GetIntruderTips}}, \code{\link{GetResultMonophyly}}, \code{\link{GetSummaryMonophyly}}, \code{\link{MonophylyData}}, \code{\link{PlotMonophyly}}
}
