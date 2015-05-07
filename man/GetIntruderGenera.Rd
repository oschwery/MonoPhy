\name{GetIntruderGenera}
\alias{GetIntruderGenera}
\alias{GetIntruderTips}
\title{
Get genera (taxa) or tips intruding a taxon from AssessMonophyly
}
\description{
Prints names of genera (or whichever taxonomic unit used) or tip names (species) interfering with the monophyly of one or several (or all) groups
}
\usage{
GetIntruderGenera(solution, genera = NULL)
GetIntruderTips(solution, genera = NULL)
}\arguments{
  \item{solution}{
Object with saved output of the 'AssessMonophyly' function.
}
  \item{genera}{
Vector containing genus names (or names of whichever taxonomic unit was used).
}
}
\details{
Can be used after 'AssessMonophyly' is run to extract intruder genera (taxa) or intruder tips from it. The argument 'genera' allows to limit the output to one or several genera (taxa) of interest.
}
\value{
Character string or list of character strings
}
\author{
Orlando Schwery
}
\note{
Can also be accessed by simply typing 'solution$Genera' or 'solution$Genera$genusname' and 'solution$Species' or 'solution$Species$genusname' respectively.
}
\seealso{
\code{\link{AssessMonophyly}}, \code{\link{GetAncNodes}}, \code{\link{MonoPhy-package}}
}
\examples{
data(Ericactree)
solution <- AssessMonophyly(Ericactree, verbosity=5)
GetIntruderGenera(solution=solution)
GetIntruderGenera(solution=solution, genera=c("Vaccinium"))
GetIntruderTips(solution=solution, genera=c("Phyllodoce", "Vaccinium", "Erica"))
}
