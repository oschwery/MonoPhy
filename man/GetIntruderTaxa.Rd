\name{GetIntruderTaxa}
\alias{GetIntruderTaxa}
\alias{GetIntruderTips}
\alias{GetOutlierTaxa}
\alias{GetOutlierTips}
\title{
Get taxa or tips which are intruders or outliers to a taxon from AssessMonophyly
}
\description{
Prints names of taxa (genera or whichever taxonomic unit used) or tip names (species, most likely) interfering with the monophyly of one or several (or all) groups, either as intruders (sharing a clade with said taxon) or outliers (being placed outside of the taxon's core clade)
}
\usage{
GetIntruderTaxa(solution, taxa = NULL, taxlevels='ALL')
GetIntruderTips(solution, taxa = NULL, taxlevels='ALL')
GetOutlierTaxa(solution, taxlevels='ALL')
GetOutlierTips(solution, taxa = NULL, taxlevels='ALL')
}\arguments{
  \item{solution}{
Object with saved output of the 'AssessMonophyly' function.
}
  \item{taxa}{
Vector containing taxon names (genus or whichever taxonomic unit was used).
}
  \item{taxlevels}{
Either an integer corresponding to the desired taxonomic level (i.e. the number of its column in the taxonimy table, not counting the tip names), or 'ALL'.  	
}
}
\details{
Can be used after 'AssessMonophyly' is run to extract intruder/outlier taxa or intruder/outlier tips from it. The argument 'taxa' allows to limit the output to one or several taxa of interest. It is not available for GetOutlierTaxa, since this is just one vector with names anyway. The argument 'taxlevels' allows to limit the output if several taxonomic levels were used; 'ALL' is default, a single level can be selected by entering its number instead. Outliers can of course only be retrieved if the argument outliercheck was set to TRUE when running AssessMonophyly.
}
\value{
List of character strings
}
\author{
Orlando Schwery
}
\note{
Noteworthy difference: intruder tips for e.g. Vaccinium will be tips of OTHER taxa that are placed inside the Vaccinium clade, whereas outlier tips for Vaccinium will be Vaccinium species that are placed outside of the Vaccinium core clade. Can also be accessed by simply typing 'solution$Taxlevel_i$IntruderTaxa' or 'solution$Taxlevel_i$IntruderTaxa$taxonname' and 'solution$Taxlevel_i$IntruderTips' or 'solution$Taxlevel_i$IntruderTips$taxonname' respectively; same for outliers.
}
\seealso{
\code{\link{AssessMonophyly}}, \code{\link{GetAncNodes}}, \code{\link{MonoPhy-package}}
}
\examples{
data(Ericactree)
solution <- AssessMonophyly(Ericactree, verbosity=5)
GetIntruderTaxa(solution=solution)
GetIntruderTaxa(solution=solution, taxa=c("Vaccinium"))
GetIntruderTips(solution=solution, taxa=c("Phyllodoce", "Vaccinium", "Erica"))
}
