\name{AssessMonophyly}
\alias{AssessMonophyly}
\title{
Assesses monophyly (or lack of it) of taxonomic groups in a phylogeny.
}
\description{
Requires rooted resolved phylogeny as input and creates a table of genera, their monophyly-status, which taxa cause problems in monophyly etc. Output can be accessed with related functions and a plot function can be used to plot the results in a number of ways.
}
\usage{
AssessMonophyly(tree, taxonomy = NULL, verbosity = 5)
}
\arguments{
  \item{tree}{
An object of type 'phy', rooted and fully resolved phylogeny. If tip labels are in the format 'genus species' or 'genus_species', the function can extract the genus names and check their monophyly. If tip labels are in another format or if the monophyly of other taxonomic groups should be tested, a taxonomy file (see 'taxonomy') is required.
}
  \item{taxonomy}{
A data frame (e.g. an imported .csv file) with at least two columns, one row per tip in tree and no header (unless specified when loading the file). Column one contains the tip labels, column two and higher the names of the taxonomic units the respective tip belongs to. The order of tip names in the file can be different from the order of the tip labels in the tree, but they have to contain the exact same names. Default is NULL.
}
  \item{verbosity}{
An integer, default is 5. Determines how many intruder taxa should be listed by name in the result table (it will list up to this many names and then add 'and X more').
}
}
\details{
The function uses \code{\link[ape]{getMRCA}} from ape and \code{\link[phytools]{getDescendants}} from phytools to determine if a genus or other taxonomic group is monophyletic or not and subsets the two to determine which taxa cause non-monophyly. From the result object, different output items (see Values) can be accessed using a set of related functions and the result can be visualized using \code{\link{PlotMonophyly}}.
}
\value{
The output of the function is a list containing the results as lists and data frames. If several taxonomic levels are analyzed, each will get a list-level containing its respective results. The different objects contained are:
  \item{Genera }{List of the names of the taxonomic groups assessed, each containing a character string with the names of genera (or whichever taxonomic level was used) which are intruders (i.e. interfere with the monophyly of that group)}
  \item{result }{Data frame containing the main results. Rows are the taxonomic groups used, columns are 'Monophyly' (yes, no or singleton), 'MRCA' (node number of inferred ancestor), 'Delta-Tips' (number of tips which share this ancestral node but do not belong to the same group), '#Intruders'(how many other groups interfere with the monophyly of the group in quesiton) and 'Intruders' (names of interfering groups, how many of these are written out is determined by the argument 'verbosity'). The whole data frame can be accessed using GetResultMonophyly}
  \item{Species }{List of the names of the taxonomic groups assessed, each containing a character string with the names of species (tip labels) which are intruders (i.e. interfere with the monophyly of that group)}
  \item{summary }{Data frame containing summary information such as number of tips, number of monophyletic genera, etc. Can be accessed using GetSummaryMonophyly.}
  \item{TipStates }{Data frame with the columns 'Tips', 'Genus' and 'Status', containing the tip labels, their associated genus (or whichever group used) and monophyly status (monophyletic, non-monophyletic or intruder) respectively.}
}
\author{
Orlando Schwery
}
\seealso{
\code{\link{GetAncNodes}}, \code{\link{GetIntruderTaxa}}, \code{\link{GetIntruderTips}}, \code{\link{GetResultMonophyly}}, \code{\link{GetSummaryMonophyly}}, \code{\link{MonophylyData}}, \code{\link{PlotMonophyly}}, \code{\link{MonoPhy-package}}
}
\examples{
data(Ericactree)
solution <- AssessMonophyly(Ericactree, verbosity=5)
GetSummaryMonophyly(solution)

#use custom taxonomic level
data(Ericactree)
data(Ericactribes)
solutiontribes <- AssessMonophyly(Ericactree, taxonomy=Ericactribes, verbosity=5)
GetSummaryMonophyly(solutiontribes)
}
