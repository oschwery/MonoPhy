\name{PlotMonophyly}
\alias{PlotMonophyly}
\title{
Plot output of AssessMonophyly
}
\description{
Allows plotting different visualisations of the results of AssessMonophyly, based on its output
}
\usage{
PlotMonophyly(solution, tree, taxlevels=1, plot.type='monophyly',
monocoll=FALSE, ladderize=TRUE, PDF=FALSE, PDF_filename='Monophylyplot.pdf',
PDF_width='auto', PDF_height='auto', mono.colour='PRGn',
tax.colour='rainbow', intrud.colour='rainbow', edge.width=3, cex=0.2,
adj.names='auto', adj.tips=0.5, label.offset='auto', type='phylogram', ...)
}
\arguments{
  \item{solution}{
Object with saved output of the 'AssessMonophyly' function.
}
  \item{tree}{
An object of type 'phy', rooted and fully resolved phylogeny (the same which was used to obtain the solution before).
}
  \item{taxlevels}{
An integer corresponding to the desired taxonomic level (i.e. the number of its column in the taxonimy table, not counting the tip names). Default is 1.  	
}
  \item{plot.type}{
Either "monophyly", "taxonomy", "intruders" or "monoVStax", see details. Default is "monophyly".
}
  \item{monocoll}{
  	If TRUE will collapse all monophyletic taxa to one single tip labeled with the taxon name. Default is FALSE.
}
  \item{ladderize}{
If TRUE will ladderize tree before reconstruction and plotting. Default is TRUE.
}
  \item{PDF}{
If TRUE will print a pdf with lenght adjusted to number of tips instead of plotting in R. Default is FALSE.
}
  \item{PDF_filename}{
Used to customize the name of pdf file created when PDF=TRUE. Default is 'Monophylyplot.pdf'.
}
  \item{PDF_width}{
Used to customize the width of pdf file created when PDF=TRUE. Default is 'auto', which will scale the width automatically based on the size of the tree.
}
  \item{PDF_height}{
Used to customize the height of pdf file created when PDF=TRUE. Default is 'auto', which will scale the height automatically based on the size of the tree.
}
  \item{mono.colour}{
Colour vector for plot type "monophyly" (and left side of "monoVStax"). Default is "PRGn", which colours monophyletic taxa green, non-monophyletic (intruded) taxa light purple and intruders purple using colours as specified in the corresponding ColorBrewer palette. Other predefined palettes are 'RdBu' (blue and red), 'PuOr' (purple and orange), 'PiYG' (green and pink) and 'BrBG' (petrol and brown). Alternatively, the colours can be customized as a vector in the order: c('monophyletic', 'non-monophyletic/intruded', 'intruder/outlier').
}
  \item{tax.colour}{
Colour vector for plot type "taxonomy" (and right side of "monoVStax"). Default is "rainbow", which assigns every taxon a different colour using the command 'rainbow()'.
}
  \item{intrud.colour}{
Colour vector for plot type "intruder". Default is "rainbow", which assigns every intruder taxon a different colour using the command 'rainbow()', while monophyletic taxa are gray and non-monophyletic (intruded) taxa are black.
}
  \item{edge.width}{
Used to customize the thickness of the tree's branches. Default is 3.
}
  \item{cex}{
Used to customize the font size of the tip labels. Default is 0.2.
}
  \item{adj.names}{
Used to change the justification of the tip label text (0: left, 0.5: centered, 1: right). Default is 'auto', which will left-justify the names in all plot types except for 'monoVStax', in which they will be centered.
}
  \item{adj.tips}{
Used to change the justification of the tip label symbols (0: left, 0.5: centered, 1: right). Default is 0.5 (centered).
}
  \item{label.offset}{
Used to adjust the distance between the tips of the tree and the tip label text. Default is 'auto', which will scale the distance based on the size of the tree.
}
  \item{type}{
Used to determine the type of phylogeny, can either be "phylogram", "cladogram", "fan" or "radial". The type "unrooted" is not available, since the assessment of monophyly requires a rooted tree. The default is "phylogram".
}
  \item{...}{
Other plot.phylo arguments.
}
}
\details{
Using the output object of \code{\link{AssessMonophyly}}, the result can be visualized in a couple of ways. Specifying the 'plot.type' "monophyly" plots a tree colouring the branches according to whether the respective clades are monophyletic, non-monophyletic or intruders/outliers. 'taxonomy' simply colour codes the different taxa. 'intruders' colours monophyletic groups gray, non-monophyletic groups black and the intruders and outliers according to the taxonomic group they belong to. Finally, 'monoVStax', plots two trees in a mirrored fashion, with the 'monophyly' type on the left and the 'taxonomy' type on the right.
It is important to note that this is merely supposed to be a useful visualization and not a biologically meaningful reconstruction.
The colour vectors must be of sufficient length if customized, i.e. three colours for mono.colour, as many as taxa for tax.colour and as many as intruder taxa for intrud.colour.
Collapsing monophyletic taxa by specifying 'monocoll = TRUE' allows to focus on the problem zones and should especially be useful for larger phylogenies.
}
\note{
Currently, if outliers are checked for (by specifying outliercheck as TRUE when running 'AssessMonophyly'), the plotting function will treat outliers and intruders the same way. This may be changed in future versions.
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
solution <- AssessMonophyly(Ericactree)

PlotMonophyly(solution=solution, tree=Ericactree, plot.type='monophyly', ladderize=TRUE,
mono.colour='PRGn')
PlotMonophyly(solution=solution, tree=Ericactree, plot.type='taxonomy', ladderize=TRUE)

# especially for larger phylogenies it is recommended to print the plot to a pdf file instead of
#plotting, for easier inspection
PlotMonophyly(solution=solution, tree=Ericactree, plot.type='monoVStax', ladderize=TRUE,
PDF=TRUE, mono.colour='PRGn', PDF_filename='MonoTaxplot.pdf')
PlotMonophyly(solution=solution, tree=Ericactree, plot.type='taxonomy', ladderize=TRUE,
PDF=TRUE, PDF_filename='Taxplot.pdf')
}
