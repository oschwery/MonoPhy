# Code to collapse monophyletic taxa into one branch (exized from plotting function)
# written by Orlando Schwery 2. May 2022
MonoColl <-
function(solution, tree, taxlevels=1, ladderize=TRUE) {
# warnings and data preparation
    if (taxlevels == 'ALL') {  # test if taxlevels argumetn has correct format and display error if not
		stop(" 'ALL' is not an option for plotting!")
    }
    if (inherits(taxlevels, 'numeric') & taxlevels > length(solution)) {  # check if specified taxlevel is among the ones in the solution object and display error if not
		stop('Requested taxonomic level not available (less levels specified as analysis input)!')
    }

	if (ladderize == TRUE) {  # ladderizes the tree before starting, if specified
        tree <- ladderize(tree)
    }
    tip.states <- solution[[taxlevels]]$TipStates  # extract tip states data frame of requested taxlevel from solution
    tip.states <- tip.states[match(tree$tip.label, tip.states$Tip), ]  # match/sort order of tip states in data frame with order of tip labels in tree
	row.names(tip.states) <- 1:nrow(tip.states)  # renumber rows of sorted table
# collapsing monophyletic groups to one tip per group
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
	tip.states.temp <- tip.states  # create temporary copy of tip.states table
	for (icollstates in 1:length(colltipse)) {  # loop through vector of tips to discard...
	    tip.states.temp <- tip.states.temp[!(tip.states.temp$Tip == tree$tip.label[colltipse[icollstates]]), ]  # and drop them from temporary table
	}
	row.names(tip.states.temp) <- 1:nrow(tip.states.temp)  # renumber rownames of new table
	for (inameadj in 1:length(tip.states.temp$Tip)) {  #loop through tip names of table...
		if (solution[[taxlevels]]$result[tip.states.temp$Taxon[inameadj], "Monophyly"] == "Yes") {  # ...and check if they belong to a monophyletic group (according to results table in solution)...
			tip.states.temp$Tip[inameadj] <- tip.states.temp$Taxon[inameadj]  # ...and then replace tip name with taxon name
		}
	}
	newlabels <- tree$tip.label  # create temporary vector of tip labels
	for (itips in 1:length(keeptips)) {  # loop through vector of tips to keep
	    newlabels[keeptips[itips]] <- as.character(tip.states[keeptips[itips], "Taxon"])  # replace tip names of tips to keep with the name of their taxon
	}
	tree$tip.label <- newlabels  # replace tip labels of tree with modified ones
	tree <- drop.tip(tree, colltipse)  # drop discard tips from tree
	tip.states <- tip.states.temp  # replace tip states table with modified one
    return(list(phy=tree, tip.states=tip.states))
}
