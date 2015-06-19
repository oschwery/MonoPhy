# assesses monophyly of genera (or customized units) and makes the result available in different ways (tables, onjects, plot...)
# written by Orlando Schwery 2015

AssessMonophyly <-
function(tree, taxonomy=NULL, verbosity=5, outliercheck=TRUE, outlierlevel=0.5, taxizelevel= NULL, taxizedb='both', taxizepref='ncbi') {
# initial tests and data preparation
if (!is.binary.tree(tree)) {  # checks and returns error if tree is not bifurcating
    stop('Phylogeny is not strictly bifurcating/resolved!')
}
if (!is.rooted(tree)) {  # checks and returns error if tree is not rooted
    stop('Phylogeny must be rooted!')
}

if (is.null(taxonomy)){  # extract list of genera from tree's tip labels
    for (i in 1:length(tree$tip.label)) {
        if (grepl(("_| "), tree$tip.label[i]) == FALSE){  # checks if genus and species epithet of tip labels are separated by space or underscore and returns error if not
            stop('Tip labels do not contain underscore separating genus name from species epithet!')
        }
    }
    
    f <- function(s) strsplit(s, ("_| "))[[1]][1]  # function with split criteria: split names at underscore and keep first part (genus name)
    split.taxa <- sapply(tree$tip.label, f)  # apply split functon to tree
    taxa <- as.vector(unique(split.taxa))  # create vector of genera in tree without duplicates
    taxsetnames <- c('taxa')
    } else {
        if (!is.null(taxonomy) && taxonomy != 'taxize') {  # use loaded taxonomy file
            if (length(taxonomy[, 1]) != length(tree$tip.label)) {  # checks and returns error if taxonomy file has more entries than tree has tips
                stop('Number of rows of taxonomy file is not equal to number of taxa (note: table should not have a header)!')
            }
            if (length(taxonomy[1, ]) < 2) {  # checks and returns error if taxonomy file doesn't have at least two columns
                stop('Taxonomy file needs at least 2 columns: tip labels and taxonomic group!')
            }
            taxchecktree <- c()
            for (itaxcheck in 1:length(tree$tip.label)) {
                taxchecktree <- c(taxchecktree, tree$tip.label[itaxcheck] %in% taxonomy[, 1])
            }
            taxintruderstree <- c()
            if ('FALSE' %in% taxchecktree) {
                positionstree <- grep('FALSE', taxchecktree)
                taxintruderstree <- c(taxintruderstree, tree$tip.label[positionstree])
                message(paste('\n'), appendLF=TRUE)
                message(paste('Tip-labels which do not occur in taxonfile:', '[', length(taxintruderstree), '/', length(tree$tip.label), ']', collapse=" "), appendLF=TRUE)
                message(paste(taxintruderstree, collapse=", "), appendLF=TRUE)
                message(paste('\n'), appendLF=TRUE)
            }
            taxcheckfile <- c()
            for (itaxcheck2 in 1:length(taxonomy[, 1])) {
                taxcheckfile <- c(taxcheckfile, taxonomy[itaxcheck2, 1] %in% tree$tip.label)
            }
            taxintrudersfile <- c()
            if ('FALSE' %in% taxcheckfile) {
                positionsfile <- grep('FALSE', taxcheckfile)
                taxintrudersfile <- c(taxintrudersfile,as.character(taxonomy[positionsfile, 1]))
                message(paste('Taxon names in file which do not occur in tip-labels of tree:', '[', length(taxintrudersfile), '/', length(taxonomy[, 1]), ']', collapse=" "), appendLF=TRUE)
                message(paste(taxintrudersfile, collapse=", "), appendLF=TRUE)
                message(paste('\n'), appendLF=TRUE)
            }
            if ('FALSE' %in% (taxchecktree)) {
                stop('The taxon names of tree and taxonfile do not match (see above)!')
            }
            if ('FALSE' %in% (taxcheckfile)) {
                stop('The taxon names of tree and taxonfile do not match (see above)!')
            }
            taxsetnames <- c()
            taxsets <- list()
            for (jtax in 1:(length(taxonomy[1, ])-1)){
                nametax <- paste("taxa", jtax, sep = "")
                taxsetnames <- c(taxsetnames, nametax)
                tmp <- as.vector(unique(taxonomy[, (jtax)+1]))  # if all is correct, makes vector taxonomic units (without doubles) 
                taxsets[[nametax]] <- tmp
            }
        } else {
            if (taxonomy == 'taxize') {  # build taxonomy file from web ressources using taxize
                taxafromweb <- tax_name(tree$tip.label, get=taxizelevel, db=taxizedb, pref=taxizepref)  # get taxonomy data
                taxafromwebtable <- matrix(data = NA, nrow=length(tree$tip.label),ncol=ncol(taxafromweb)+1)  #build empty matrix
                taxafromwebtable[, 1] <- tree$tip.label  # add tip names from tree
                for (iweb in 1:ncol(taxafromweb)) {  #add acquired taxon names for tips for each taxonomic level acquired
                    taxafromwebtable[, iweb+1] <- taxafromweb[, iweb]
                }
                taxafromwebtable[is.na(taxafromwebtable)] <- "unknown"
                taxonomy <- as.data.frame(taxafromwebtable)  # turn matrix into data frame and feed to further function
            }  
        }
    }
# actual assessment
finallist <- list()
for (ifullround in 1:length(taxsetnames)){  # Assess monophyly for every taxon set used
    if (is.null(taxonomy)){  # assing correct taxonset
        taxa <- taxa
    } else {
        taxa <- unlist(taxsets[ifullround])
        taxa <- taxa[!taxa %in% c("unknown", NA)]
    }
    # create empty objects to be filled by function
    intruder.genus <- list()  # list of genera causing non-monophyly
    intruder.genus.full <- c()  # vector of ALL general causing non-monophyly
    intruder.species <- list()  # list of species causing non-monophyly
    intruder.species.full <- c()  # vector of ALL species causing non-monophyly
    intruder.names <- c()  # names for intruder sub-lists
    if (outliercheck == TRUE) {
        outlist.summary <- matrix(NA, nrow=6, ncol=3)  # final output summary matrix
        dfheaders <- c("Taxon", "Monophyly", "MRCA", "#Tips", "Delta-Tips", "#Intruders", "Intruders", "#Outliers", "Outliers")  # headers for 'outlist'
        outlist <- matrix(NA, nrow=length(taxa), ncol=9)  # final output matrix
        outlier.genus <- list()  # list of genera causing non-monophyly as outliers
        outlier.genus.full <- c()  # vector of ALL general causing non-monophyly as outliers
        outlier.species <- list()  # list of species causing non-monophyly as outliers
        outlier.species.full <- c()  # vector of ALL species causing non-monophyly as outliers
        outlier.names <- c()  # names for outlier sub-lists
    } else {
        outlist.summary <- matrix(NA, nrow=5, ncol=3)  # final output summary matrix
        dfheaders <- c("Taxon", "Monophyly", "MRCA", "#Tips", "Delta-Tips", "#Intruders", "Intruders")  # headers for 'outlist'
        outlist <- matrix(NA, nrow=length(taxa), ncol=7)  # final output matrix
    }
    tip.states.matrix <- matrix(NA, nrow=length(tree$tip.label), ncol=3)  # states for plotting
    
    # loop assessing monophyly
    for (i in 1:length(taxa)) {  # for every genus in the tree
        if (is.null(taxonomy)){  # genera extracted from tip labels if no taxonomy file loaded
            ancnode <- getMRCA(tree, tip=c(tree$tip.label[c(grep(taxa[i], tree$tip.label))]))  # determine Most Recent Common Ancestor for all taxa of genus
        } else {  # units taken from file if loaded
            subtips <- subset(taxonomy, as.character(taxonomy[, (ifullround+1)]) == as.character(taxa[i]))  # get tips associated with group
            subtipsnr <- c()
            for (sbts in 1: nrow(subtips)) {
                sbtname <- subtips[sbts,1]
                sbtnr <- which(tree$tip.label==sbtname)
                subtipsnr <- c(subtipsnr, sbtnr)
            }
            ancnode <- getMRCA(tree, tip=c(subtipsnr))  # get MRCA for tips associated with group
        }
        if (length(ancnode) == 0) { # if monotypic i.e. only tip of given group
            if (outliercheck == TRUE) {
                outlist[i, ] <- c(taxa[i], "Monotypic", "NA", 1, "NA", "NA", "", "NA", "")  # UPDATE OUTPUT MATRIX, mark as monotypic if only one tip for this genus
            } else {
                outlist[i, ] <- c(taxa[i], "Monotypic", "NA", 1, "NA", "NA", "")  # UPDATE OUTPUT MATRIX, mark as monotypic if only one tip for this genus
            }
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
                    if (outliercheck == TRUE) {
                        outlist[i, ] <- c(taxa[i], "Yes", ancnode, length(taxtips), "0", "0", "", "NA", "")  # UPDATE OUTPUT MATRIX, mark as monophyletic
                    } else {
                        outlist[i, ] <- c(taxa[i], "Yes", ancnode, length(taxtips), "0", "0", "")  # UPDATE OUTPUT MATRIX, mark as monophyletic
                    }
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
                                    
                                    nodechoice <- which(c(length(taxtips1), length(taxtips2))==max(c(length(taxtips1), length(taxtips2))))
                                    if (length(nodechoice) > 1) {
                                        nodechoice2 <- which(c((length(taxtips1)/length(anctips1)), (length(taxtips2)/length(anctips2)))==max(c((length(taxtips1)/length(anctips1)), (length(taxtips2)/length(anctips2)))))
                                        if (nodechoice2 == 1) {
                                            subtaxtips <- taxtips1
                                            subancnames <- ancnames1
                                            start.node <- daughter1
                                        }
                                        if (nodechoice2 == 2) {
                                            subtaxtips <- taxtips2
                                            subancnames <- ancnames2
                                            start.node <- daughter2
                                        }
                                        if (length(nodechoice) > 1) {
                                            subtaxtips <- c(taxtips1, taxtips2)
                                            subancnames <- c(ancnames1, ancnames2)
                                            start.node <- parent.node
                                            break
                                        }
                                    }
                                    if (nodechoice == 1) {
                                        subtaxtips <- taxtips1
                                        subancnames <- ancnames1
                                        start.node <- daughter1
                                    }
                                    if (nodechoice == 2) {
                                        subtaxtips <- taxtips2
                                        subancnames <- ancnames2
                                        start.node <- daughter2
                                    }
                                    tiplevels <- length(subtaxtips)/length(subancnames)
                                }

                                outlier.tips <- setdiff(taxtips, subtaxtips)
                                
                                if (length(outlier.tips) != 0) {
                                    outlier.species <- c(outlier.species, list(Tips=outlier.tips))  # update list of intruder tip labels
                                    outlier.species.full <- c(outlier.species.full, outlier.tips)  # update vector of ALL intruder species
                                    outlier.names <- c(outlier.names, taxa[i])  # update names vector for intruders
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
    taxcount.monoph <- sapply(levels(outframe$Monophyly),function(x){sum(as.numeric(outframe$'#Tips'[outframe$Monophyly==x]))},USE.NAMES=F)
    taxcount.frame <- data.frame(status=levels(outframe$Monophyly), taxcount.monoph)
    rownames(taxcount.frame) <- taxcount.frame[, 1]
    taxcount.frame[, 1] <- NULL
    if (outliercheck == TRUE) {
        outlist.summary[, 2] <- c(length(taxa), countframe["Yes","Freq"], countframe["No","Freq"], countframe["Monotypic","Freq"], length(intruder.genus.all), length(outlier.names))  # populate summary table with counts of respective groups
        outlist.summary[, 3] <- c(length(tree$tip.label), taxcount.frame["Yes", "taxcount.monoph"], taxcount.frame["No", "taxcount.monoph"], countframe["Monotypic","Freq"], length(intruder.species.all), length(outlier.species.all))
    } else {
        outlist.summary[, 2] <- c(length(taxa), countframe["Yes","Freq"], countframe["No","Freq"], countframe["Monotypic","Freq"], length(intruder.genus.all))  # populate summary table with counts of respective groups
        outlist.summary[, 3] <- c(length(tree$tip.label), taxcount.frame["Yes", "taxcount.monoph"], taxcount.frame["No", "taxcount.monoph"], countframe["Monotypic","Freq"], length(intruder.species.all))
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
finallist
}