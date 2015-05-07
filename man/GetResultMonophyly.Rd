\name{GetResultMonophyly}
\alias{GetResultMonophyly}
\alias{GetSummaryMonophyly}
\title{
Get main result table from 'AssessMonophyly' output
}
\description{
Prints data frame with main result table or containing summary information such as number of tips, number of monophyletic genera, etc. from the output of 'AssessMonophyly'.
}
\usage{
GetResultMonophyly(solution)
GetSummaryMonophyly(solution)
}
\arguments{
  \item{solution}{
Object with saved output of the 'AssessMonophyly' function.
}
}
\details{
Can be used after 'AssessMonophyly' is run to extract the results or summary table from it. In the results table, rows are the taxonomic groups used, columns are 'Monophyly' (yes, no or singleton), 'MRCA' (node number of inferred ancestor), 'Delta-Tips' (number of tips which share this ancestral node but do not belong to the same group), '#Intruders'(how many other groups interfere with the monophyly of the group in quesiton) and 'Intruders' (names of interfering groups, how many of these are written out is determined by the argument 'verbosity').
The summary table contains number of genera (taxa), tips (species), monophyletic and non-monophyletic groups, singletons, intruder genera and intruder tips.
}
\value{
Data frame containing the main results/summary.
}
\author{
Orlando Schwery
}
\note{
Can also be accessed by simply typing 'solution$result' or 'solution$summary' respectively.
}
\seealso{
\code{\link{AssessMonophyly}}, \code{\link{MonoPhy-package}}
}
\examples{
data(Ericactree)
solution <- AssessMonophyly(Ericactree, verbosity=5)
GetSummaryMonophyly(solution)
GetResultMonophyly(solution)
}
