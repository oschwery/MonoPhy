\name{GetAncNodes}
\alias{GetAncNodes}
\title{
Get MRCA nodes of taxonomic groups from 'AssessMonophyly' output
}
\description{
Print MRCA node numbers of taxonomic groups from the output of 'AssessMonophyly', either the whole list or for specific groups.
}
\usage{
GetAncNodes(solution, genera = NULL)
}
\arguments{
  \item{solution}{
Object with saved output of the 'AssessMonophyly' function.
}
  \item{genera}{
Vector containing genus names (or names of whichever taxonomic unit was used).
}
}
\details{
Can be used after 'AssessMonophyly' is run to extract MRCA nodes of groups from it. The argument 'genera' allows to limit the output to one or several genera (taxa) of interest.
}
\value{
Data frame.
}
\author{
Orlando Schwery
}
\seealso{
\code{\link{AssessMonophyly}}}, \code{\link{MonoPhy-package}}
\examples{
data(Ericactree)
solution <- AssessMonophyly(Ericactree, verbosity=5)
GetAncNodes(solution=solution, genera=c("Phyllodoce", "Vaccinium", "Erica"))
}
