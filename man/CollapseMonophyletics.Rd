\name{CollapseMonophyletics}
\alias{CollapseMonophyletics}
\title{
Get tree with collapsed monophyletic groups
}
\description{
Creates a new tree object with all monophyletic groups collapsed to one tip each, based on the output of 'AssessMonophyly'.
}
\usage{
CollapseMonophyletics(solution, tree, taxlevels = 1, ladderize = TRUE)
}
\arguments{
  \item{solution}{
Object with saved output of the 'AssessMonophyly' function.
}
  \item{tree}{
An object of type 'phy', rooted phylogeny (the same which was used to obtain the solution before).
}
  \item{taxlevels}{
An integer or name corresponding to the desired taxonomic level (i.e. the number of its column in the taxonomy table, not counting the tip names). Default is 1.
}
  \item{ladderize}{
If TRUE will ladderize tree before collapsing clades. Default is TRUE.
}
}
\details{
Can be used after 'AssessMonophyly' is run to obtain a tree object where all groups that were inferred as monophyletic at a given taxonomic level will be collapsed to single tips and relabeled according to the name of their group. The resulting tree will have the same topology as seen when plotting the original tree with setting 'monocoll=TRUE' in 'PlotMonophyly'.
}
\value{
A tree object of type 'phy'.
}
\author{
Orlando Schwery, Peter Cowman
}
\note{
The code is largely the same as 'PlotMonophyly' uses to create a plot with collapsed monophyletic groups.
}
\seealso{
\code{\link{AssessMonophyly}}, \code{\link{PlotMonophyly}}, \code{\link{MonoPhy-package}}
}
\examples{
data(Ericactree)
solution <- AssessMonophyly(Ericactree, verbosity=5)
collapsetree <- CollapseMonophyletics(solution, Ericactree, taxlevels = 1, ladderize = TRUE)
}
