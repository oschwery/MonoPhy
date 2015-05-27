\name{MonophylyData}
\docType{data}
\alias{MonophylyData}
\alias{Ericactree}
\alias{Ericactribes}
\alias{Ericacsubfams}
\title{
Example dataset for the package MonoPhy.
}
\description{
Example dataset for the package MonoPhy. One phylogeny and two custom taxonomy files.
}
\usage{
data(Ericactree)
data(Ericactribes)
data(Ericacsubfams)
}
\details{
Ericactree is a phylogeny of the angiosperm family Ericaceae, as presented in Schwery et al. (2014), pruned to 77 species. Ericactribes is an example file for how to specify which taxonomic groups should be tested for monophyly (instead of genus names taken from the tip labels of the tree), here exemplified by tribes. In Ericacsubfams we have two columns, one with tribes and one with subfamilies. For the sake of the example, some taxa have been intentionally mislabelled in the taxonomy files: the tribes of the tips 'Vaccinium_vitis_ideae' and 'Zenobia_pulverulenta' have been swapped and 'Rhodothamnus_chamaecystus' was labelled to be a member of the 'Vaccinioideae'. 
}
\references{
Schwery, O., Onstein, R.E., Bouchenak-Khelladi, Y., Xing, Y., Carter, R.J. and Linder, H.P. (2014), As old as the mountains: the radiations of the Ericaceae. New Phytologist. doi: 10.1111/nph.13234
}
\seealso{
\code{\link{MonoPhy-package}}
}
\keyword{datasets}
