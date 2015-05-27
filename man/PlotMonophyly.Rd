\name{PlotMonophyly}
\alias{PlotMonophyly}
\title{
Plot output of AssessMonophyly
}
\description{
Allows plotting different visualisations of the results of AssessMonophyly, based on its output
}
\usage{
function(solution, tree, taxlevels=1, type='monophyly', ladderize=TRUE, PDF=FALSE, PDF_filename='Monophylyplot.pdf', mono.colour='PRGn', tax.colour='rainbow', intrud.colour='rainbow', edge.width=3, cex=0.2, ...)
}
\arguments{
  \item{solution}{
Object with saved output of the 'AssessMonophyly' function.
}
  \item{tree}{
An object of type 'phy', rooted and fully resolved phylogeny.
}
  \item{taxlevels}{
An integer corresponding to the desired taxonomic level (i.e. the number of its column in the taxonimy table, not counting the tip names). Default is 1.  	
}
  \item{type}{
Either "monophyly", "taxonomy", "intruders" or "monoVStax", see details. Default is "monophyly".
}
  \item{ladderize}{
If TRUE will ladderize tree before reconstruction and plotting. Default is TRUE.
}
  \item{PDF}{
If TRUE will print a pdf with lenght adjusted to number of tips instead of plotting in R. Default is FALSE.
}
  \item{PDF_filename}{
Used to customize the name of pdf file created when PDF=TRUE. Default is Monophylyplot.pdf.
}
  \item{mono.colour}{
Colour vector for plot type "monophyly" (and left side of "monoVStax"). Default is "PRGn", which colours monophyletic taxa green, non-monophyletic (intruded) taxa light purple and intruders purple using colours as specified in the corresponding ColorBrewer palette. Other predefined palettes are 'RdBu' (blue and red), 'PuOr' (purple and orange), 'PiYG' (green and pink) and 'BrBG' (petrol and brown). Alternatively, the colours can be customized as a vector in the order: c('monophyletic', 'non-monophyletic/intruded', 'intruder').
}
  \item{tax.colour}{
Colour vector for plot type "taxonomy" (and right side of "monoVStax"). Default is "rainbow", which assigns every taxon a different colour using the command 'rainbow()'.
}
  \item{intrud.colour}{
Colour vector for plot type "intruder". Default is "rainbow", which assigns every intruder taxon a different colour using the command 'rainbow()', while monophyletic taxa are gray and non-monophyletic (intruded) taxa are black.
}
  \item{...}{
Other plot.phylo arguments.
}
}
\details{
Using the output object of \code{\link{AssessMonophyly}}, the result can be visualized in a couple of ways. Specifying the 'type' "monophyly" plots a tree colouring the branches according to whether the respective clades are monophyletic, non-monophyletic or intruders. 'taxonomy' simply colour codes the different taxa. 'intruders' colours monophyletic groups gray, non-monophyletic groups black and the intruders according to the taxonomic group they belong to. Finally, 'monoVStax', plots two trees in a mirrored fashion, with the 'monophyly' type on the left and the 'taxonomy' type on the right.
It is important to note that this is merely supposed to be a useful visualization and not a biologically meaningful reconstruction.
The colour vectors must be of sufficient length if customized, i.e. three colours for mono.colour, as many as taxa for tax.colour and as many as intruder taxa for intrud.colour.
}
\value{
Plots phylogeny or prints it to pdf.
}
\author{
Orlando Schwery
}
\seealso{
\code{\link{AssessMonophyly}}, \code{\link{MonoPhy-package}}
}
\examples{
data(Ericactree)
solution <- AssessMonophyly(Ericactree, verbosity=5)

PlotMonophyly(solution=solution, tree=Ericactree, type='monophyly', ladderize=TRUE, mono.colour='PRGn')
PlotMonophyly(solution=solution, tree=Ericactree, type='taxonomy', ladderize=TRUE)

# especially for larger phylogenies it is recommended to print the plot to a pdf file instead of plotting, for easier inspection
PlotMonophyly(solution=solution, tree=Ericactree, type='monoVStax', ladderize=TRUE, PDF=TRUE, mono.colour='PRGn', PDF_filename='MonoTaxplot.pdf')
PlotMonophyly(solution=solution, tree=Ericactree, type='taxonomy', ladderize=TRUE, PDF=TRUE, PDF_filename='Taxplot.pdf')
}
