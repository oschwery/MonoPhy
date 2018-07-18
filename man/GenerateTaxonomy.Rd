\name{GenerateTaxonomy }
\alias{GenerateTaxonomy}
\title{
Generate a taxonomy table for all the taxa in a tree
}
\description{
Requires rooted phylogeny as input and creates a data.frame of the taxonomic hierarchy (one row per species)
}
\usage{
GenerateTaxonomy(tree, taxizepref='ncbi')
}
\arguments{
  \item{tree}{
An object of type 'phy', a rooted phylogeny. Multifurcating trees are accepted, but will be dealt with in a conservative manner (i.e. if different taxa share a multifurcation, they will be considered non-monophyletic) If tip labels are in the format 'genus_species', the function can extract the genus names and check their monophyly. If tip labels are in another format or if the monophyly of other taxonomic groups should be tested, a taxonomy file (see 'taxonomy') is required.
}
  \item{taxizepref}{
Your choice of database: could be one of 'ncbi, itis, eol, col, tropicos, gbif, nbn, worms, natserv, bold,' or 'wiki'. See ?taxize::classification for more info.
}
}
\details{
The function uses 'taxize::classification' to extract information for each species. In the case of multiple names for a species, it uses the first.
}
\value{
The output object of the function is a data.frame with one row per species on the tree and one column per rank (some species may have NA for some ranks). Note that it will return an object with many more levels than AssessMonophyly will use, so subset.
}
\author{
Brian O'Meara
}
\seealso{
\code{\link{GetAncNodes}}, \code{\link{GetIntruderTaxa}}, \code{\link{GetIntruderTips}}, \code{\link{GetOutlierTaxa}}, \code{\link{GetOutlierTips}},\code{\link{GetResultMonophyly}}, \code{\link{GetSummaryMonophyly}}, \code{\link{MonophylyData}}, \code{\link{PlotMonophyly}}, \code{\link{MonoPhy-package}}
}
