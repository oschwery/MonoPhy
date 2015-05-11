\name{GetResultMonophyly}
\alias{GetResultMonophyly}
\alias{GetSummaryMonophyly}
\title{
Get main result table from 'AssessMonophyly' output
}
\description{
Extracts data frame with main result table or containing summary information such as number of tips, number of monophyletic genera, etc. from the output of 'AssessMonophyly'.
}
\usage{
GetResultMonophyly(solution, taxlevels='ALL')
GetSummaryMonophyly(solution, taxlevels='ALL')
}
\arguments{
  \item{solution}{
Object with saved output of the 'AssessMonophyly' function.
}
  \item{taxlevels}{
Either an integer corresponding to the desired taxonomic level (i.e. the number of its column in the taxonimy table, not counting the tip names), or 'ALL'.  	
}
}
\details{
Can be used after 'AssessMonophyly' is run to extract the results or summary table from it. In the results table, rows are the taxonomic groups used, columns are 'Monophyly' (yes, no or singleton), 'MRCA' (node number of inferred ancestor), 'Delta-Tips' (number of tips which share this ancestral node but do not belong to the same group), '#Intruders'(how many other groups interfere with the monophyly of the group in quesiton) and 'Intruders' (names of interfering groups, how many of these are written out is determined by the argument 'verbosity').
The summary table contains number of genera (taxa), tips (species), monophyletic and non-monophyletic groups, singletons, intruder genera and intruder tips.
The argument 'taxlevels' allows to limit the output if several taxonomic levels were used; 'ALL' is default, a specific level can be selected by entering its corresponding number instead.
}
\value{
List of data frame(s) containing the main results/summary.
}
\author{
Orlando Schwery
}
\note{
Can also be accessed by simply typing 'solution$Taxlevel_i$result' or 'solution$Taxlevel_i$summary' respectively.
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
