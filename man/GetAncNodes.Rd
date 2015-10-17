\name{GetAncNodes}
\alias{GetAncNodes}
\title{
Get MRCA nodes of taxonomic groups from 'AssessMonophyly' output
}
\description{
Print MRCA node numbers of taxonomic groups from the output of 'AssessMonophyly', either the whole list or for specific groups.
}
\usage{
GetAncNodes(solution, taxa = NULL, taxlevels='ALL')
}
\arguments{
  \item{solution}{
Object with saved output of the 'AssessMonophyly' function.
}
  \item{taxa}{
Vector containing taxon names (genera or whichever taxonomic unit was used). Default is 'NULL'
}
  \item{taxlevels}{
Either an integer corresponding to the desired taxonomic level (i.e. the number of its column in the taxonimy table, not counting the tip names), or 'ALL' (which is the default).  	
}
}
\details{
Can be used after 'AssessMonophyly' is run to extract MRCA nodes of taxa from it. The argument 'taxa' allows to limit the output to one or several taxa of interest. The argument 'taxlevels' allows to limit the output if several taxonomic levels were used; 'ALL' is default, a specific level can be selected by entering its number instead.
}
\value{
List of data frames.
}
\author{
Orlando Schwery
}
\seealso{
\code{\link{AssessMonophyly}}, \code{\link{MonoPhy-package}}
}
\examples{
data(Ericactree)
solution <- AssessMonophyly(Ericactree)
GetAncNodes(solution=solution, taxa=c("Phyllodoce", "Vaccinium", "Erica"))
}
