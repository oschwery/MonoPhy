# get new tree object with monophyletic groups collapsed to one tip each, based on the output of AssessMonophyly
# written by Orlando Schwery and Peter Cowman 2016

CollapseMonophyletics <- function(solution, tree, taxlevels = 1, ladderize = TRUE) {  # collapsing monophyletic groups to one tip per group
if (ladderize == TRUE) {  # ladderizes the tree before starting, if specified
      tree <- ladderize(tree)
  }
  tip.states <- solution[[taxlevels]]$TipStates  # extract tip states data frame of requested taxlevel from solution
  tip.states <- tip.states[match(tree$tip.label, tip.states$Tip), ]  # match/sort order of tip states in data frame with order of tip labels in tree
  row.names(tip.states) <- 1:nrow(tip.states)  # renumber rows of sorted table
  tip.states$Tip <- as.character(tip.states$Tip)  # turn tip labels into characters
  tip.states$Taxon <- as.character(tip.states$Taxon)  # turn taxon names into characers
  tip.states$Status <- as.character(tip.states$Status)  # turn status into characters
  alltaxa <- as.vector(unique(tip.states[, "Taxon"]))  # create vector of all unique taxon names in table
  keeptips <- c()  # create empty vector for tips to keep
  colltipse <- c()  # create empty vector for tips to collapse/discard
  for (icoll in 1:length(alltaxa)) {  # loop through unique taxa
      if (solution[[taxlevels]]$result[alltaxa[icoll], "Monophyly"] == "Yes") {  # if the current taxon is monophyletic (according to results table in solution)
      matchtips <- which(tip.states[, "Taxon"] == alltaxa[icoll])  # get vector with all tip.states rows that belong (match) to current taxon
      colltips <- c(matchtips[2:length(matchtips)])  # select all but first of these tips...
      colltipse <- c(colltipse, colltips)  # ...and add them to discard vector
      keeptip <- c(matchtips[1])  # select first of these tips...
      keeptips <- c(keeptips, keeptip)  # ...and add it to vector with tips to keep
      }
  }
  newlabels <- tree$tip.label  # create temporary vector of tip labels
  for (itips in 1:length(keeptips)) {  # loop through vector of tips to keep
      newlabels[keeptips[itips]] <- as.character(tip.states[keeptips[itips], "Taxon"])  # replace tip names of tips to keep with the name of their taxon
  }
  tree$tip.label <- newlabels  # replace tip labels of tree with modified ones
  coltree <- drop.tip(tree, colltipse)  # drop discard tips from tree
}
