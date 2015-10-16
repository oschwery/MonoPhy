# assesses monophyly of genera (or customized units) and makes the result available in different ways (tables, onjects, plot...)
# written by Orlando Schwery 2015
AssessMonophyly <-
function(tree, taxonomy=NULL, verbosity=5, outliercheck=TRUE, outlierlevel=0.5, taxizelevel= NULL, taxizedb='both', taxizepref='ncbi', taxask=FALSE, taxverbose=FALSE) {
# initial tests and data preparation
if (!is.binary.tree(tree)) {  # checks and returns error if tree is not bifurcating
    stop('Phylogeny is not strictly bifurcating/resolved!')
}
if (!is.rooted(tree)) {  # checks and returns error if tree is not rooted
    stop('Phylogeny must be rooted!')
}
if (is.null(taxonomy)){  # if no taxonomy file is specified, extract list of genera from tree's tip labels
    for (i in 1:length(tree$tip.label)) {  # loop through tip labels of tree
        if (grepl(("_| "), tree$tip.label[i]) == FALSE){  # checks if genus and species epithet of tip labels are separated by space or underscore and returns error if not
            stop('Tip labels do not contain underscore separating genus name from species epithet!')
        }
    }
    f <- function(s) strsplit(s, ("_| "))[[1]][1]  # function with split criteria: split names at underscore and keep first part (genus name)
    split.taxa <- sapply(tree$tip.label, f)  # apply split functon to tree
    taxa <- as.vector(unique(split.taxa))  # create vector of genera in tree without duplicates
    taxsetnames <- c('taxa')  # assign 'taxa' as taxsetnames
    } else {
        if (!is.null(taxonomy) && taxonomy != 'taxize') {  # if argument 'taxonomy' is not NULL and not taxize, use loaded taxonomy file
            if (length(taxonomy[, 1]) != length(tree$tip.label)) {  # checks and returns error if taxonomy file has more or less entries than tree has tips
                stop('Number of rows of taxonomy file is not equal to number of taxa (note: table should not have a header)!')
            }
            if (length(taxonomy[1, ]) < 2) {  # checks and returns error if taxonomy file doesn't have at least two columns
                stop('Taxonomy file needs at least 2 columns: tip labels and taxonomic group!')
            }
            taxchecktree <- c()  # create empty vector to be filled with presence of tip labels in taxfile
            for (itaxcheck in 1:length(tree$tip.label)) {  # loop through tip labels
                taxchecktree <- c(taxchecktree, tree$tip.label[itaxcheck] %in% taxonomy[, 1])  # check for every name in tree if it is present in taxonomy file
            }
            taxintruderstree <- c()  # create empty vector to be filled with tip labels that are abesent in taxfile
            if ('FALSE' %in% taxchecktree) {  # if there are any missing names in taxonomy file...
                positionstree <- grep('FALSE', taxchecktree)  # ...get their position...
                taxintruderstree <- c(taxintruderstree, tree$tip.label[positionstree])  # ...add tip labels of that position to vector
                message(paste('\n'), appendLF=TRUE)
                message(paste('Tip-labels which do not occur in taxonomy file:', '[', length(taxintruderstree), '/', length(tree$tip.label), ']', collapse=" "), appendLF=TRUE)  # display numbers of missing tip labels in taxonomy file
                message(paste(taxintruderstree, collapse=", "), appendLF=TRUE)  # display names of missing tips
                message(paste('\n'), appendLF=TRUE)
            }
            taxcheckfile <- c()  # create empty vector to be filled with presence of taxfile names in tip labels
            for (itaxcheck2 in 1:length(taxonomy[, 1])) {  # loop through names column of taxfile
                taxcheckfile <- c(taxcheckfile, taxonomy[itaxcheck2, 1] %in% tree$tip.label)  # check for every name in taxonomy file if it is present in tip labels of the tree
            }
            taxintrudersfile <- c()  # create empty vector to be filled with taxfile names that are abesent in tip labels
            if ('FALSE' %in% taxcheckfile) {  # if there are any missing names in tip labels...
                positionsfile <- grep('FALSE', taxcheckfile)  # ...get their position...
                taxintrudersfile <- c(taxintrudersfile,as.character(taxonomy[positionsfile, 1]))  # ...add names of that position to vector
                message(paste('Taxon names in file which do not occur in tip-labels of tree:', '[', length(taxintrudersfile), '/', length(taxonomy[, 1]), ']', collapse=" "), appendLF=TRUE)  # display numbers of missing names in tip labels
                message(paste(taxintrudersfile, collapse=", "), appendLF=TRUE)  # display names of missing taxa
                message(paste('\n'), appendLF=TRUE)
            }
            if ('FALSE' %in% (taxchecktree)) {  # if missing names in file, stop and display error
                stop('The taxon names of tree and taxonfile do not match (see above)!')
            }
            if ('FALSE' %in% (taxcheckfile)) {  # if missing names in tree, stop and display error
                stop('The taxon names of tree and taxonfile do not match (see above)!')
            }
            taxsetnames <- c()  # create empty vector to fill with names for taxsets
            taxsets <- list()  # create empty list to fill with taxonomic units
            for (jtax in 1:(length(taxonomy[1, ]) - 1)) {  # loop through taxon file
                nametax <- paste("taxa", jtax, sep = "")  # create taxon name label
                taxsetnames <- c(taxsetnames, nametax)  # add label to names vector
                tmp <- as.vector(unique(taxonomy[, (jtax)+1]))  # if all is correct, makes vector of taxonomic units (without doubles) 
                taxsets[[nametax]] <- tmp  # add names vector to taxets list, labelled with taxon name label
            }
        } else {  # if not NULL and not taxfile but taxize
            if (taxonomy == 'taxize') {  # build taxonomy file from web ressources using taxize
                taxafromweb <- tax_name(tree$tip.label, get=taxizelevel, db=taxizedb, pref=taxizepref, ask=taxask, verbose=taxverbose)  # get taxonomy data from web
                taxafromwebtable <- matrix(data = NA, nrow=length(tree$tip.label),ncol=ncol(taxafromweb)-1)  #build empty matrix with dimensions by number of tips and retrieved taxonomic data
                taxafromwebtable[, 1] <- tree$tip.label  # add tip names from tree
                if (length(unique(taxafromweb[, 3])) == 1 & is.na(unique(taxafromweb[, 3])[1])) {  # check if there was any information retrieved and display error if not
                    stop('There was no data found for any of the tips!')
                }
                if (nrow(taxafromweb) > length(tree$tip.label)) {  # check if more record entires retrieved than tips in tree
                    taxafromweb <- unique(taxafromweb[, 2:3])  # if entires from two databases are the same, delete one
                    rownames(taxafromweb) <- c(1:nrow(taxafromweb))  # renumber rows
                }
                if (nrow(taxafromweb) > length(tree$tip.label)) {  # check again if more record entires retrieved than tips in tree
                    temp <- taxafromweb  # copy table to temporary object
                    droppers <- c()  # create empty vector to be filled with entires to be dropped
                    counterweb <- table(temp[, 1], useNA='ifany')  # create counting table with number of entries per taxon in retrieved data table
                    for (iwebtax in 1:nrow(temp)) {  # loop through retrieved table
                        if (counterweb[temp[iwebtax, 1]] > 1) {  # check if more than one entry per species is retrieved (according to count table)
                            if (is.na(temp[iwebtax, 2])) {  # if the duplicate currently being looked at did not retrieve a result...
                                droppers <- c(droppers, iwebtax)  # ...add it to droplist
                            }
                        }
                        temp2 <- temp[-droppers, ]  # remove drop entires
                        rownames(temp2) <- c(1:nrow(temp2))  # renumber rows
                        taxafromweb <- temp2  # replace retrieved table with cleaned up version
                    }
                }
                for (iweb in 1:(ncol(taxafromweb) - 2)) {  #add acquired taxon names for tips for each taxonomic level acquired
                    taxafromwebtable[, iweb + 1] <- taxafromweb[, iweb + 1]  # add retrieved entries to matrix
                    rownames(taxafromweb) <- c(1:nrow(taxafromweb))  # renumber rownames
                }
                taxafromwebtable[is.na(taxafromwebtable)] <- "unknown"  # replace all NAs with "unknown"
                taxonomy <- as.data.frame(taxafromwebtable)  # turn matrix into data frame and feed to further function
                taxsetnames <- c()  # create empty vector to be filled with taxonomy set names
                taxsets <- list()  # create empty list to be filled with taxonomy lists 
                for (jtax in 1:(length(taxonomy[1, ]) - 1)) {  # loop through taxonomy file
                    nametax <- paste("taxa", jtax, sep = "")  # create taxon name label
                    taxsetnames <- c(taxsetnames, nametax)  # add label to names vector
                    tmp <- as.vector(unique(taxonomy[, (jtax)+1]))  # if all is correct, makes vector taxonomic units (without doubles) 
                    taxsets[[nametax]] <- tmp  # add names vector to taxets list, labelled with taxon name label
                }
            }  
        }
    }
# actual assessment
finallist <- list()  # create empty list to be filled with final results
for (ifullround in 1:length(taxsetnames)){  # Assess monophyly for every taxon set used
    if (is.null(taxonomy)){  # if no taxonomy file specified...
        taxa <- taxa  # ...assign taxon list
    } else {  # if taxonomy file or taxize...
        taxa <- unlist(taxsets[ifullround])  # ...assign taxon list of current taxlevel...
        taxa <- taxa[!taxa %in% c("unknown", NA)]  # ...and remove NA's and unknown taxa
    }
# create empty objects to be filled by function
    intruder.genus <- list()  # empty list for genera causing non-monophyly
    intruder.genus.full <- c()  # empty vector for ALL general causing non-monophyly
    intruder.species <- list()  # empty list for species causing non-monophyly
    intruder.species.full <- c()  # empty vector for ALL species causing non-monophyly
    intruder.names <- c()  # empty vector for names for intruder sub-lists
    if (outliercheck == TRUE) {  # if outliers should be checked for (add additional objects and columns/rows
        outlist.summary <- matrix(NA, nrow=6, ncol=3)  # create final output summary matrix
        dfheaders <- c("Taxon", "Monophyly", "MRCA", "#Tips", "Delta-Tips", "#Intruders", "Intruders", "#Outliers", "Outliers")  # headers for 'outlist'
        outlist <- matrix(NA, nrow=length(taxa), ncol=9)  # final output matrix
        outlier.genus <- list()  # list of genera causing non-monophyly as outliers
        outlier.genus.full <- c()  # vector of ALL general causing non-monophyly as outliers
        outlier.species <- list()  # list of species causing non-monophyly as outliers
        outlier.species.full <- c()  # vector of ALL species causing non-monophyly as outliers
        outlier.names <- c()  # names for outlier sub-lists
    } else {  # if no outliers being checked for
        outlist.summary <- matrix(NA, nrow=5, ncol=3)  # final output summary matrix
        dfheaders <- c("Taxon", "Monophyly", "MRCA", "#Tips", "Delta-Tips", "#Intruders", "Intruders")  # headers for 'outlist'
        outlist <- matrix(NA, nrow=length(taxa), ncol=7)  # final output matrix
    }
    tip.states.matrix <- matrix(NA, nrow=length(tree$tip.label), ncol=3)  # states for plotting
# loop assessing monophyly
    for (i in 1:length(taxa)) {  # loop through every genus in the tree
        if (is.null(taxonomy)) {  # genera extracted from tip labels if no taxonomy file loaded
            ancnode <- getMRCA(tree, tip=c(tree$tip.label[c(grep(paste("^", taxa[i], "_", sep=""), tree$tip.label))]))  # determine Most Recent Common Ancestor for all taxa of genus
        } else {  # units taken from file if loaded
            subtips <- subset(taxonomy, as.character(taxonomy[, (ifullround+1)]) == as.character(taxa[i]))  # get tips associated with current group
            subtipsnr <- c()  # create empty vector to be filled with tip numbers
            for (sbts in 1: nrow(subtips)) {  # looping through all subtips assigned to current group
                sbtname <- subtips[sbts, 1]  # extract name
                sbtnr <- which(tree$tip.label == sbtname)  # extract number of subtip
                subtipsnr <- c(subtipsnr, sbtnr)  # add subtip nr to numbers vector
            }
            ancnode <- getMRCA(tree, tip=c(subtipsnr))  # get MRCA for tips associated with group
        }
        if (length(ancnode) == 0) { # if monotypic i.e. only tip of given group
            if (outliercheck == TRUE) {  # if outliers are checked for (more columns)
                outlist[i, ] <- c(taxa[i], "Monotypic", "NA", 1, "NA", "NA", "", "NA", "")  # UPDATE OUTPUT MATRIX, mark as monotypic if only one tip for this genus
            } else {  # if outliers are not checked for (less columns)
                outlist[i, ] <- c(taxa[i], "Monotypic", "NA", 1, "NA", "NA", "")  # UPDATE OUTPUT MATRIX, mark as monotypic if only one tip for this genus
            }
        } else {  # if not monotypic
            anctips <- getDescendants(tree, ancnode)  # determine all descendants of previously determined MRCA
            ancnames <- tree$tip.label[c(anctips)]  # extract names of those descendants
            ancnames <- ancnames[!is.na(ancnames)]  # ommit NA's (caused by descendants which are internal nodes and not tips)
            if (is.null(taxonomy)) {  # genera extracted from tip labels if no taxonomy file loaded
                taxtips <- tree$tip.label[c(grep(paste("^",taxa[i],"_", sep=""), tree$tip.label))]  # get tip names of genus in question
            } else {  # if taxonomy file loaded
                taxtips <- subtips[, 1]  # get vector of tip names of genus in question
            }
            if (length(ancnames) == length(taxtips)) {  # determine if all MRCA descendants = genus members. Genus is monophyletic if yes.
                if (outliercheck == TRUE) {  # if outliers are checked for (more columns)
                    outlist[i, ] <- c(taxa[i], "Yes", ancnode, length(taxtips), "0", "0", "", "NA", "")  # UPDATE OUTPUT MATRIX, mark as monophyletic
                } else {  # if outliers are checked for (less columns)
                    outlist[i, ] <- c(taxa[i], "Yes", ancnode, length(taxtips), "0", "0", "")  # UPDATE OUTPUT MATRIX, mark as monophyletic
                }
            } else {  # if taxon is not monophyletic
                intruder.tips <- setdiff(ancnames, taxtips)  # determine intruders tip labels, i.e. descendants of MRCA which are not genus members
                if (is.null(taxonomy)){  # get intruder genus names if no taxonomy file loaded
                            f2 <- function(s) strsplit(s, ("_| "))[[1]][1]  # function with split criteria: split names at underscore and keep first part (genus name)
                            split.taxa2 <- sapply(intruder.tips, f2)  # apply split function to intruder tip labels
                            intruder.taxa <- as.vector(unique(split.taxa2))  # create vector of intruder genera
                        } else {  # get intruder taxonomic levels from taxonomy file
                            subtaxa <- c()
                            for (j in 1:length(intruder.tips)) {
                                subtaxon <- rbind(subset(taxonomy, taxonomy[, 1] == intruder.tips[j]))  # extract taxon for each intruder tip...
                                subtaxa <- rbind(subtaxa, subtaxon)  # ... and add them up
                            }                      
                            intruder.taxa <- as.vector(unique(subtaxa[, ifullround+1]))  # create vector of intruder taxa
                        }
                        outlier.tips <- c()
                        if (outliercheck == TRUE) {  # distinguish outliers if TRUE
                            tiplevels <- length(taxtips)/length(ancnames)
                            if (tiplevels < outlierlevel ) {  # check if meeting criteria
                                start.node <- ancnode
                                while (tiplevels < outlierlevel) {  # search for subclade that meets criteria
                                    subtaxtips <- c()  # reset taxtips
                                    subancnames <- c()  # reset ancnames
                                    parent.node <- start.node # set parent node
                                    daughter.nodes <- Children(tree, parent.node) # find direct descendant nodes
                                    daughter1 <- daughter.nodes[1]
                                    daughter2 <- daughter.nodes[2]
                                    
                                    anctips1 <- getDescendants(tree, daughter1)  # determine all descendants of daughter1
                                    ancnames1 <- tree$tip.label[c(anctips1)]  # extract names of those descendants
                                    ancnames1 <- ancnames1[!is.na(ancnames1)]  # ommit NA's (caused by descendants which are internal nodes and not tips)
                                    taxtips1 <- intersect(taxtips, ancnames1)  # get taxon members of subclade1
                                    
                                    anctips2 <- getDescendants(tree, daughter2)  # determine all descendants of daughter2
                                    ancnames2 <- tree$tip.label[c(anctips2)]  # extract names of those descendants
                                    ancnames2 <- ancnames2[!is.na(ancnames2)]  # ommit NA's (caused by descendants which are internal nodes and not tips)
                                    taxtips2 <- intersect(taxtips, ancnames2)  # get taxon members of subclade2
                                    # pick from the daughter nodes
                                    nodechoice <- which(c(length(taxtips1), length(taxtips2))==max(c(length(taxtips1), length(taxtips2))))  # determine daughter with more tips of focal taxon
                                    if (length(nodechoice) > 1) {  # if equal number of taxon tips in both daughers:
                                        nodechoice2 <- which(c((length(taxtips1)/length(anctips1)), (length(taxtips2)/length(anctips2)))==max(c((length(taxtips1)/length(anctips1)), (length(taxtips2)/length(anctips2)))))  # determine daugther with higher ratio of focal taxon
                                        if (length(nodechoice2) > 1) {  # if equal ratio of taxon tips in both daughters: keep both daughters as core clade
                                            subtaxtips <- c(taxtips1, taxtips2)
                                            subancnames <- c(ancnames1, ancnames2)
                                            start.node <- parent.node
                                            break
                                        } else if (nodechoice2 == 1) {  # if daughter 1 chosen: set as new start-node and -clade
                                            subtaxtips <- taxtips1
                                            subancnames <- ancnames1
                                            start.node <- daughter1
                                        } else if (nodechoice2 == 2) {  # if daughter 2 chosen: set as new start-node and -clade
                                            subtaxtips <- taxtips2
                                            subancnames <- ancnames2
                                            start.node <- daughter2
                                        }
                                    } else if (nodechoice == 1) {  # if daughter 1 chosen: set as new start-node and -clade
                                        subtaxtips <- taxtips1
                                        subancnames <- ancnames1
                                        start.node <- daughter1
                                    } else if (nodechoice == 2) {  # if daughter 2 chosen: set as new start-node and -clade
                                        subtaxtips <- taxtips2
                                        subancnames <- ancnames2
                                        start.node <- daughter2
                                    }
                                    tiplevels <- length(subtaxtips)/length(subancnames)  # reassess status of current clade
                                }
                                if (tiplevels < 1) {  # if intruders are present, check if early-diverging
                                    EDtaxtips1 <- c()
                                    EDtaxtips2 <- c()
                                    #while (length(EDtaxtips1) == 0 | length(EDtaxtips2) == 0) {  # search for node whose daughers both include members of the focal taxon
                                    repeat{
                                        EDparent.node <- start.node # set parent node
                                        if (EDparent.node <= length(tree$tip.label)) {
                                            subtaxtips <- tree$tip.label[EDparent.node]
                                            subancnames <-tree$tip.label[EDparent.node]
                                            start.node <- EDparent.node
                                            break
                                        }
                                        EDdaughter.nodes <- Children(tree, EDparent.node)  # find direct descendant nodes
                                        EDdaughter1 <- EDdaughter.nodes[1]
                                        EDdaughter2 <- EDdaughter.nodes[2]
                                        
                                        EDanctips1 <- getDescendants(tree, EDdaughter1)  # determine all descendants of EDdaughter1
                                        EDancnames1 <- tree$tip.label[c(EDanctips1)]  # extract names of those descendants
                                        EDancnames1 <- EDancnames1[!is.na(EDancnames1)]  # ommit NA's (caused by descendants which are internal nodes and not tips)
                                        EDtaxtips1 <- intersect(taxtips, EDancnames1)  # get taxon members of subclade1
                                        
                                        EDanctips2 <- getDescendants(tree, EDdaughter2)  # determine all descendants of EDdaughter2
                                        EDancnames2 <- tree$tip.label[c(EDanctips2)]  # extract names of those descendants
                                        EDancnames2 <- EDancnames2[!is.na(EDancnames2)]  # ommit NA's (caused by descendants which are internal nodes and not tips)
                                        EDtaxtips2 <- intersect(taxtips, EDancnames2)  # get taxon members of subclade2
                                        # select node to continue
                                        if (length(EDtaxtips1) == 0) {  # if EDdaugher1 has no decendant of focal taxon, continue with EDdaughter2
                                            start.node <- EDdaughter2
                                        } else if (length(EDtaxtips2) == 0) {  # if EDdaugher2 has no decendant of focal taxon, continue with EDdaughter1
                                            start.node <- EDdaughter1
                                        } else {
                                            subtaxtips <- c(EDtaxtips1, EDtaxtips2)
                                            subancnames <- c(EDancnames1, EDancnames2)
                                            start.node <- EDparent.node
                                            break
                                        }
                                    }
                                }
                                outlier.tips <- setdiff(taxtips, subtaxtips)  # determine outliers
                                
                                if (length(outlier.tips) != 0) {
                                    outlier.species <- c(outlier.species, list(Tips=outlier.tips))  # update list of outlier tip labels
                                    outlier.species.full <- c(outlier.species.full, outlier.tips)  # update vector of ALL outlier species
                                    outlier.names <- c(outlier.names, taxa[i])  # update names vector for outliers
                                }
                                
                                intruder.tips <- setdiff(subancnames, subtaxtips)
                                
                                if (is.null(taxonomy)){  # get intruder genus names if no taxonomy file loaded
                                     f2 <- function(s) strsplit(s, ("_| "))[[1]][1]  # function with split criteria: split names at underscore and keep first part (genus name)
                                    split.taxa2 <- sapply(intruder.tips, f2)  # apply split function to intruder tip labels
                                    intruder.taxa <- as.vector(unique(split.taxa2))  # create vector of intruder genera
                                } else {  # get intruder taxonomic levels from taxonomy file
                                    subtaxa <- c()
                                    for (j in 1:length(intruder.tips)) {
                                        subtaxon <- rbind(subset(taxonomy, taxonomy[, 1] == intruder.tips[j]))  # extract taxon for each intruder tip...
                                        subtaxa <- rbind(subtaxa, subtaxon)  # ... and add them up
                                    }                      
                                    intruder.taxa <- as.vector(unique(subtaxa[, ifullround+1]))  # create vector of intruder taxa
                                }
                            }
                        }
                        if (length(intruder.taxa) != 0) {
                            intruder.genus <- c(intruder.genus, list(Taxa=intruder.taxa))  # update list of intruder genera
                            intruder.genus.full <- c(intruder.genus.full,intruder.taxa)  # update vector of ALL intruder genera
                            intruder.species <- c(intruder.species, list(Tips=intruder.tips))  # update list of intruder tip labels
                            intruder.species.full <- c(intruder.species.full,intruder.tips)  # update vector of ALL intruder species
                            intruder.names <- c(intruder.names, taxa[i])  # update names vector for intruders
                        }
                        if (outliercheck == TRUE) {
                            if (length(intruder.taxa) <= verbosity) {  # give verbose intruder list if less than 5 intruding genera
                                if (length(outlier.tips) <= verbosity) {
                                    outlist[i, ] <- c(taxa[i], "No", ancnode, length(taxtips), (length(ancnames) - length(taxtips)), length(intruder.taxa), paste(intruder.taxa, collapse=", "), length(outlier.tips), paste(outlier.tips, collapse=", "))  # UPDATE OUTPUT MATRIX, mark as non-monophyletic and list intruder genera
                                } else {
                                   outlist[i, ] <- c(taxa[i], "No", ancnode, length(taxtips), (length(ancnames) - length(taxtips)), length(intruder.taxa), paste(intruder.taxa, collapse=", "), length(outlier.tips), paste(outlier.tips[1], "and", (length(outlier.tips) - 1), "more.", collapse=", "))  # UPDATE OUTPUT MATRIX, mark as non-monophyletic and list intruder genera
                                }
                            } else {
                                if (length(outlier.tips) <= verbosity) {
                                    outlist[i, ] <- c(taxa[i], "No", ancnode, length(taxtips), (length(ancnames) - length(taxtips)), length(intruder.taxa), paste(intruder.taxa[1], "and", (length(intruder.taxa) - 1), "more.", collapse=", "), length(outlier.tips), paste(outlier.tips, collapse=", "))  # UPDATE OUTPUT MATRIX, mark as non-monophyletic and list intruder genera
                                } else {
                                    outlist[i, ] <- c(taxa[i], "No", ancnode, length(taxtips), (length(ancnames) - length(taxtips)), length(intruder.taxa), paste(intruder.taxa[1], "and", (length(intruder.taxa) - 1), "more.", collapse=", "), length(outlier.tips), paste(outlier.tips[1], "and", (length(outlier.tips) - 1), "more.", collapse=", "))  # UPDATE OUTPUT MATRIX, mark as non-monophyletic and list intruder genera
                                }
                            }
                        } else {
                            if (length(intruder.taxa) <= verbosity) {  # give verbose intruder list if less than 5 intruding genera
                                outlist[i, ] <- c(taxa[i], "No", ancnode, length(taxtips), (length(ancnames) - length(taxtips)), length(intruder.taxa), paste(intruder.taxa, collapse=", "))  # UPDATE OUTPUT MATRIX, mark as non-monophyletic and list intruder genera
                            } else {
                                outlist[i, ] <- c(taxa[i], "No", ancnode, length(taxtips), (length(ancnames) - length(taxtips)), length(intruder.taxa), paste(intruder.taxa[1], "and", (length(intruder.taxa) - 1), "more.", collapse=", "))  # UPDATE OUTPUT MATRIX, mark as non-monophyletic and list intruder genera
                            }
                        }
                }
            }
        }
    # prepare outputs
    intruder.genus.all <- unique(intruder.genus.full)  # vector of ALL intruder genera without doubles
    intruder.species.all <- unique(intruder.species.full)  # vector of ALL intruder species without doubles
    outlier.species.all <- c()
    if (outliercheck == TRUE) {
        outlier.species.all <- unique(outlier.species.full)  # vector of ALL outlier species without doubles
    }
    outframe <- data.frame(outlist)  # turn final output matrix into data frame
    names(outframe) <- dfheaders  # apply headers
    rownames(outframe) <- outframe[, 1]  # assign first column as row names
    outframe[, 1] <- NULL  # delete first column (since now row names)
    names(intruder.genus) <- intruder.names  # apply names (of intruded genera) to intruder genus lists
    names(intruder.species) <- intruder.names  # apply names (of intruded genera) to intruder species lists
    
    if (is.null(taxonomy)){
        tip.states.matrix[, 1] <- tree$tip.label  # adding species names to matrix
        f3 <- function(s) strsplit(s, ("_| "))[[1]][1]  # function with split criteria: split names at underscore and keep first part (genus name)
        tip.states.matrix[, 2] <- sapply(tip.states.matrix[, 1], f3)  # apply split function, create genus name column
    } else {
        tip.states.matrix[, 1] <- as.vector(taxonomy[, 1])  # adding species names to matrix
        tip.states.matrix[, 2] <- as.vector(taxonomy[, ifullround+1])  # create taxon name column
    }
    for (i in 1:length(tree$tip.label)){  #for every matrix entry
        if (tip.states.matrix[i, 2] == "unknown" | is.na(tip.states.matrix[i, 2])) {
            tip.states.matrix[i, 3] <- "unknown"
        } else {
            if (tip.states.matrix[i, 1] %in% intruder.species.all == TRUE){  # score species as intruder if in global intruder list
                tip.states.matrix[i, 3] <- "Intruder"
            } else {
                if (tip.states.matrix[i, 1] %in% outlier.species.all == TRUE){  # score species as intruder if in global intruder list
                    tip.states.matrix[i, 3] <- "Outlier"
                } else {
                    if (outframe[tip.states.matrix[i, 2], "Monophyly"] == "Monotypic"){  # if not intruder, score as monophyletic if monotypic
                        tip.states.matrix[i, 3] <- "Monophyletic"
                    } else {
                        if (outframe[tip.states.matrix[i, 2], "Monophyly"] == "Yes"){  # if not intruder nor monotypic, score as monophyletic if monophyletic
                            tip.states.matrix[i, 3] <- "Monophyletic"
                        } else {
                            if (outframe[tip.states.matrix[i, 2], "Monophyly"] == "No"){  # if not intruder nor monotypic, score as monophyletic if monophyletic
                                tip.states.matrix[i, 3] <- "Non-Monophyletic"  # if not intruder, monotypic or monophyletic, score as non-monophyletic
                            }
                        }
                    }      
                }
            }
        } 
    }
    tip.states.frame <- as.data.frame(tip.states.matrix)  # turn into data frame
    #rownames(tip.states.frame) <- tip.states.matrix[, 1]  # assign tip labels as row names
    #tip.states.frame[, 1] <- NULL  # delete first row, since now row names
    colnames(tip.states.frame) <- c("Tip", "Taxon", "Status")  # assign column names
    if (outliercheck == TRUE) {
        names(outlier.species) <- outlier.names  # apply names (of intruded genera) to intruder species lists
        outlist.summary[, 1] <- c("Total", "Monophyletic", "Non-Monophyletic", "Monotypic", "Intruder", "Outlier")  # row names for summary table
    } else {
        outlist.summary[, 1] <- c("Total", "Monophyletic", "Non-Monophyletic", "Monotypic", "Intruder")  # row names for summary table
    }
    counttable <- table(outframe[, "Monophyly"])  # tabulate monophyly results
    countframe <- as.data.frame(counttable)  # turn into data frame
    rownames(countframe) <- countframe[, 1]  # assign first column as row names
    countframe[, 1] <- NULL  # delete first column (since now row names)

    mono.count <- c()
    nonmono.count <-c()

    for (itaxcount1 in 1:length(taxa)) {
        if (outlist[itaxcount1,2] == "Yes") {
            mono.count <- c(mono.count, as.numeric(outlist[itaxcount1,4]))
        }
    }
    for (itaxcount2 in 1:length(taxa)) {
        if (outlist[itaxcount2,2] == "No") {
            nonmono.count <- c(nonmono.count, as.numeric(outlist[itaxcount2,4]))
        }
    }
    if (outliercheck == TRUE) {
        outlist.summary[, 2] <- c(length(taxa), countframe["Yes","Freq"], countframe["No","Freq"], countframe["Monotypic","Freq"], length(intruder.genus.all), length(outlier.names))  # populate summary table with counts of respective groups
        outlist.summary[, 3] <- c(length(tree$tip.label), sum(mono.count), sum(nonmono.count), countframe["Monotypic","Freq"], length(intruder.species.all), length(outlier.species.all))
    } else {
        outlist.summary[, 2] <- c(length(taxa), countframe["Yes","Freq"], countframe["No","Freq"], countframe["Monotypic","Freq"], length(intruder.genus.all))  # populate summary table with counts of respective groups
        outlist.summary[, 3] <- c(length(tree$tip.label), sum(mono.count), sum(nonmono.count), countframe["Monotypic","Freq"], length(intruder.species.all))
    }
    outframe.summary <- data.frame(outlist.summary)  # turn final output matrix summary into data frame
    rownames(outframe.summary) <- outframe.summary[, 1]  # assign first column as row names
    outframe.summary[, 1] <- NULL  # delete first column (since now row names)
    colnames(outframe.summary) <- c("Taxa", "Tips")
    if (outliercheck == TRUE) {
        outputlist <- list(IntruderTaxa=intruder.genus, IntruderTips=intruder.species, OutlierTaxa=outlier.names, OutlierTips=outlier.species, result=outframe, summary=outframe.summary, TipStates=tip.states.frame) # concatenate intruder lists, final output table and summmary to one list object
    } else {
        outputlist <- list(IntruderTaxa=intruder.genus, IntruderTips=intruder.species, result=outframe, summary=outframe.summary, TipStates=tip.states.frame) # concatenate intruder lists, final output table and summmary to one list object
    }
    nameout <- paste("Taxlevel", ifullround, sep = "_")  # name for output subsection
    finallist[[nameout]] <- outputlist #add list for this round of the loop to final list
    }
    finallist  # return final outputlist
}