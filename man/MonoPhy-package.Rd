\name{MonoPhy-package}
\alias{MonoPhy-package}
\alias{MonoPhy}
\docType{package}
\title{
Allows to explore monophyly (or lack of it) of taxonomic groups in a phylogeny.
}
\description{
Requires rooted phylogeny as input and creates a table of taxa, their monophyly-status, which other taxa cause problems in monophyly (as intruders or outliers) etc. Different information can be extracted from the output and a plot function allows visualization of the results in a number of ways.
}
\details{
\tabular{ll}{
Package: \tab MonoPhy\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2016-03-08\cr
License: \tab GPL-3\cr
}
The package allows to explore monophyly in phylogenetic trees in a quick and simple  way. AssessMonophyly determines whether taxa in a phylogeny are monophyletic. GetSummaryMonophyly extracts a summary table from the results, GetResultMonophyly extracts the main result table, GetIntruderTaxa, GetIntruderTips, GetOutlierTaxa and GetOutlierTips extract names of taxa which interfere with the monophyly of another or the same taxon respectively and GetAncNodes a table of MRCA nodes. PlotMonophyly finally, allows visualizing the results in a couple of different ways.
}
\author{
Orlando Schwery

Maintainer: Orlando Schwery <oschwery@vols.utk.edu>
}
\seealso{
\code{\link{AssessMonophyly}}, \code{\link{GetAncNodes}}, \code{\link{GetIntruderTaxa}}, \code{\link{GetIntruderTips}}, \code{\link{GetOutlierTaxa}}, \code{\link{GetOutlierTips}} \code{\link{GetResultMonophyly}}, \code{\link{GetSummaryMonophyly}}, \code{\link{MonophylyData}}, \code{\link{PlotMonophyly}}
}
