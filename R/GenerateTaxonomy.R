GenerateTaxonomy <- function(tree, taxizepref="ncbi") {
  species <- tree$tip.label
  taxonomies <- taxize::classification(species, db=taxizepref, rows=1)
  bad.ones <- c()
  for (i in sequence(length(taxonomies))) {
    if (class(taxonomies[[i]]) == "logical") {
      bad.ones <- c(bad.ones,i)
    } else {
      tax.vector <- taxonomies[[i]]$name
      names(tax.vector) <- make.names(taxonomies[[i]]$rank, unique=TRUE)
      taxonomies[[i]] <- data.frame(t(rev(tax.vector)), stringsAsFactors=FALSE)
    }
  }
  if(length(bad.ones)>0) {
    warning(paste("The following taxa had no matches:", paste(bad.ones, collapse=", ")))
    taxonomies <- taxonomies[-bad.ones]
  }
  final.taxonomy <- plyr::rbind.fill(taxonomies)
  final.taxonomy[(1+nrow(final.taxonomy)):(nrow(final.taxonomy)+length(bad.ones)),1] <- gsub("_", " ", names(taxonomies)[bad.ones])
  return(final.taxonomy)
}
