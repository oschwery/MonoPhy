# assesses monophyly of genera (or customized units) and makes the result available in different ways (tables, onjects, plot...)
# written by Orlando Schwery 2015

AssessMonophyly <-
function(tree, taxonomy=NULL, verbosity=5) {

if (!is.binary.tree(tree)) {  # checks and returns error if tree is not bifurcating
    stop('phylogeny is not strictly bifurcating/resolved')
}
if (!is.rooted(tree)) {  # checks and returns error if tree is not rooted
    stop('phylogeny must be rooted')
}

if (is.null(taxonomy)){  # extract list of genera from tree's tip labels
    for (i in 1:length(tree$tip.label)) {
        if (grepl(("_| "), tree$tip.label[i]) == FALSE){  # checks if genus and species epithet of tip labels are separated by space or underscore and returns error if not
            stop('tip labels do not contain underscore/space separating genus name from species epithet')
        }
    }
    
    f <- function(s) strsplit(s, ("_| "))[[1]][1]  # function with split criteria: split names at underscore and keep first part (genus name)
    split.taxa <- sapply(tree$tip.label, f)  # apply split functon to tree
    taxa <- as.vector(unique(split.taxa))  # create vector of genera in tree without duplicates
    } else {  # use loaded file
        if (length(taxonomy[, 1]) != length(tree$tip.label)) {  # checks and returns error if taxonomy file has more entries than tree has tips
            stop('number of rows of taxonomy file is not equal to number of taxa (possibility: if table has header, header=TRUE must be specified when loading file)')
        }
        if (length(taxonomy[1, ]) < 2) {  # checks and returns error if taxonomy file doesn't have at least two columns
            stop('taxonomy file needs at least 2 columns: tip labels and taxonomic group')
        }
        taxa <- as.vector(unique(taxonomy[, 2]))  # if all is correct, makes vector taxonomic units (without doubles)
}

# create empty objects to be filled by function
outlist <- matrix(NA, nrow=length(taxa), ncol=6)  # final output matrix
intruder.genus <- list()  # list of genera causing non-monophyly
intruder.genus.full <- c()  # vector of ALL general causing non-monophyly
intruder.species <- list()  # list of species causing non-monophyly
intruder.species.full <- c()  # vector of ALL species causing non-monophyly
intruder.names <- c()  # names for intruder sub-lists
outlist.summary <- matrix(NA, nrow=7, ncol=2)  # final output summary matrix
dfheaders <- c("Taxon", "Monophyly", "MRCA", "Delta-Tips", "#Intruders", "Intruders")  # headers for 'outlist'

# loop assessing monophyly
for (i in 1:length(taxa)) {  # for every genus in the tree
    if (is.null(taxonomy)){  # genera extracted from tip labels if no taxonomy file loaded
        ancnode <- getMRCA(tree, tip=c(tree$tip.label[c(grep(taxa[i], tree$tip.label))]))  # determine Most Recent Common Ancestor for all taxa of genus
    } else {  # units taken from file if loaded
        subtips <- subset(taxonomy, as.character(taxonomy[, 2]) == taxa[i])  # get tips associated with group
        ancnode <- getMRCA(tree, tip=c(subtips[, 1]))  # get MRCA for tips associated with group
    }
    if (length(ancnode) == 0) { # if singleton i.e. only tip of given group
        outlist[i, ] <- c(taxa[i], "Singleton", "NA", "NA", "NA", "")  # UPDATE OUTPUT MATRIX, mark as singleton if only one tip for this genus
    } else {
            anctips <- getDescendants(tree, ancnode)  # determine all descendants of previously determined MRCA
            ancnames <- tree$tip.label[c(anctips)]  # extract names of those descendants
            ancnames <- ancnames[!is.na(ancnames)]  # ommit NA's (caused by descendants which are internal nodes and not tips)
            if (is.null(taxonomy)){  # genera extracted from tip labels if no taxonomy file loaded
                taxtips <- tree$tip.label[c(grep(taxa[i], tree$tip.label))]  # get tip names of genus in question
            } else {  # if taxonomy file loaded
                taxtips <- subtips[, 1]  # get vector of tip names of genus in question
            }
            if (length(ancnames) == length(taxtips)) {  # determine if all MRCA descendants = genus members. Genus is monophyletic if yes.
                outlist[i, ] <- c(taxa[i], "Yes", ancnode, "0", "0", "")  # UPDATE OUTPUT MATRIX, mark as monophyletic
                } else {
                    intruder.tips <- setdiff(ancnames, taxtips)  # determine intruders tip labels, i.e. descendants of MRCA which are not genus members
                    if (is.null(taxonomy)){  # get intruder genus names if no taxonomy file loaded
                        f2 <- function(s) strsplit(s, ("_| "))[[1]][1]  # function with split criteria: split names at underscore and keep first part (genus name)
                        split.taxa2 <- sapply(intruder.tips, f2)  # apply split function to intruder tip labels
                        intruder.taxa <- as.vector(unique(split.taxa2))  # create vector of intruder genera
                    } else {  # get intruder taxonomic levels from taxonomy file
                        subtaxa <- c()
                        for (j in 1:length(intruder.tips)) {
                            subtaxon <- rbind(subset(taxonomy, taxonomy[, 1] == intruder.tips[j]))  # extract taxon for each intruder tip...
                            subtaxa <- rbind(subtaxa, subtaxon)  # ... and add them to a vector
                        }                      
                        intruder.taxa <- as.vector(unique(subtaxa[, 2]))  # create vector of intruder taxa
                    }
                    
                    intruder.genus <- c(intruder.genus, list(Taxa=intruder.taxa))  # update list of intruder genera
                    intruder.genus.full <- c(intruder.genus.full,intruder.taxa)  # update vector of ALL intruder genera
                    
                    intruder.species <- c(intruder.species, list(Tips=intruder.tips))  # update list of intruder tip labels
                    intruder.species.full <- c(intruder.species.full,intruder.tips)  # update vector of ALL intruder species
                    intruder.names <- c(intruder.names, taxa[i])  # update names vector for intruders
                    if (length(intruder.taxa) <= verbosity) {  # give verbose intruder list if less than 5 intruding genera
                    outlist[i, ] <- c(taxa[i], "No", ancnode, (length(ancnames) - length(taxtips)), length(intruder.taxa), paste(intruder.taxa, collapse=", "))  # UPDATE OUTPUT MATRIX, mark as non-monophyletic and list intruder genera
                    } else {
                    outlist[i, ] <- c(taxa[i], "No", ancnode, (length(ancnames) - length(taxtips)), length(intruder.taxa), paste(intruder.taxa[1], "and", (length(intruder.taxa) - 1), "more.", collapse=", "))  # UPDATE OUTPUT MATRIX, mark as non-monophyletic, list first intruder genus and number of remaining ones
                    }
                }
        }
    }

# prepare outputs
intruder.genus.all <- unique(intruder.genus.full)  # vector of ALL intruder genera without doubles
intruder.species.all <- unique(intruder.species.full)  # vector of ALL intruder species without doubles

outframe <- data.frame(outlist)  # turn final output matrix into data frame
names(outframe) <- dfheaders  # apply headers
rownames(outframe) <- outframe[, 1]  # assign first column as row names
outframe[, 1] <- NULL  # delete first column (since now row names)
names(intruder.genus) <- intruder.names  # apply names (of intruded genera) to intruder genus lists
names(intruder.species) <- intruder.names  # apply names (of intruded genera) to intruder species lists
outlist.summary[, 1] <- c("Total Genera", "Total Tips", "Monophyletic", "Non-Monophyletic", "Singletons", "Intruder Genera", "Intruder Tips")  # row names for summary table
counttable <- table(outframe[, "Monophyly"])  # tabulate monophyly results
countframe <- as.data.frame(counttable)  # turn into data frame
rownames(countframe) <- countframe[, 1]  # assign first column as row names
countframe[, 1] <- NULL  # delete first column (since now row names)
outlist.summary[, 2] <- c(length(taxa), length(tree$tip.label), countframe["Yes","Freq"], countframe["No","Freq"], countframe["Singleton","Freq"], length(intruder.genus.all), length(intruder.species.all))  # populate summary table with counts of respective groups
outframe.summary <- data.frame(outlist.summary)  # turn final output matrix summary into data frame
rownames(outframe.summary) <- outframe.summary[, 1]  # assign first column as row names
outframe.summary[, 1] <- NULL  # delete first column (since now row names)

tip.states.matrix <- matrix(NA, nrow=length(tree$tip.label), ncol=3)  # states for plotting
tip.states.matrix[, 1] <- tree$tip.label  # adding species names to matrix
if (is.null(taxonomy)){
    f3 <- function(s) strsplit(s, ("_| "))[[1]][1]  # function with split criteria: split names at underscore and keep first part (genus name)
    tip.states.matrix[, 2] <- sapply(tip.states.matrix[, 1], f3)  # apply split function, create genus name column
} else {
    tip.states.matrix[, 2] <- as.vector(taxonomy[, 2])  # create taxon name column
}
for (i in 1:length(tree$tip.label)){  #for every matrix entry
    if (tip.states.matrix[i, 1] %in% intruder.species.all == TRUE){  # score species as intruder if in global intruder list
        tip.states.matrix[i, 3] <- "Intruder"
    } else {
        if (outframe[tip.states.matrix[i, 2], "Monophyly"] == "Singleton"){  # if not intruder, score as monophyletic if singleton
            tip.states.matrix[i, 3] <- "Monophyletic"
        } else {
        if (outframe[tip.states.matrix[i, 2], "Monophyly"] == "Yes"){  # if not intruder nor singleton, score as monophyletic if monophyletic
            tip.states.matrix[i, 3] <- "Monophyletic"
        } else {
            tip.states.matrix[i, 3] <- "Non-Monophyletic"  # if not intruder, singleton or monophyletic, score as non-monophyletic
        }      
        }
    }
}
tip.states.frame <- as.data.frame(tip.states.matrix)  # turn into data frame
#rownames(tip.states.frame) <- tip.states.matrix[, 1]  # assign tip labels as row names
#tip.states.frame[, 1] <- NULL  # delete first row, since now row names
colnames(tip.states.frame) <- c("Tip", "Genus", "Status")  # assign column names
outputlist <- list(Genera=intruder.genus, Species=intruder.species, result=outframe, summary=outframe.summary, TipStates=tip.states.frame) # concatenate intruder lists, final output table and summmary to one list object
}
