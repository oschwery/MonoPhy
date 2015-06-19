\name{GetResultMonophyly}
\alias{GetResultMonophyly}
\alias{GetSummaryMonophyly}
\title{
Get main result table from 'AssessMonophyly' output
}
\description{
Extracts data frame with main result table or containing summary information such as number of tips, number of monophyletic taxa, etc. from the output of 'AssessMonophyly'.
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
Can be used after 'AssessMonophyly' is run to extract the results or summary table from it. In the results table, rows are the taxonomic groups used, columns are 'Monophyly' ('Yes', 'No' or 'Monotypic'), 'MRCA' (node number of inferred ancestor), '#Tips' (number of tips assigned to this taxon),'Delta-Tips' (number of tips which share this ancestral node but do not belong to the same taxon), '#Intruders'(how many other taxa interfere with the monophyly of the taxon in quesiton) and 'Intruders' (names of interfering taxa, how many of these are written out is determined by the argument 'verbosity'). If the argument outliercheck was set to TRUE when running 'AssessMonophyly', the table will additionally include '#Outliers' (number of tips which are placed outside the core clade for that group) and 'Outliers' (names of outlier taxa, the number of which also depending on the argument 'verbosity').
The summary table contains two rows with the number of taxa (e.g. genera) and tips (e.g. species) for the total tree, and which were inferred to be monophyletic, non-monophyletic, monotypic, intruders and (if applicable) outliers.
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
