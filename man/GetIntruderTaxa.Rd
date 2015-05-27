\name{GetIntruderTaxa}
\alias{GetIntruderTaxa}
\alias{GetIntruderTips}
\title{
Get taxa or tips intruding a taxon from AssessMonophyly
}
\description{
Prints names of taxa (genera or whichever taxonomic unit used) or tip names (species, most likely) interfering with the monophyly of one or several (or all) groups
}
\usage{
GetIntruderTaxa(solution, taxa = NULL, taxlevels='ALL')
GetIntruderTips(solution, taxa = NULL, taxlevels='ALL')
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
Can be used after 'AssessMonophyly' is run to extract intruder taxa or intruder tips from it. The argument 'taxa' allows to limit the output to one or several taxa of interest. The argument 'taxlevels' allows to limit the output if several taxonomic levels were used; 'ALL' is default, a single level can be selected by entering its number instead.
}
\value{
List of character strings
}
\author{
Orlando Schwery
}
\note{
Can also be accessed by simply typing 'solution$Taxlevel_i$Taxa' or 'solution$Taxlevel_i$Taxa$taxonname' and 'solution$Taxlevel_i$Tips' or 'solution$Taxlevel_i$Tips$taxonname' respectively.
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
