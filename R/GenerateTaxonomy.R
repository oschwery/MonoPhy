GenerateTaxonomy <- function(tree, taxizepref="ncbi") {
  species <- tree$tip.label
  taxonomies <- taxize::classification(species, db=taxizedb, rows=1)
  for (i in sequence(length(taxonomies))) {
    tax.vector <- taxonomies[[i]]$name

    names(tax.vector) <- make.names(taxonomies[[i]]$rank, unique=TRUE)
    taxonomies[[i]] <- data.frame(t(rev(tax.vector)), stringsAsFactors=FALSE)
  }
  final.taxonomy <- plyr::rbind.fill(taxonomies)
  return(final.taxonomy)
}
