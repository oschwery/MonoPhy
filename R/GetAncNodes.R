#get list of MRCA nodes of desired genera only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetAncNodes <-
function(solution, taxa=NULL, taxlevels='ALL') {
    allnodes <- list()
    if (taxlevels=='ALL') {
        for (i in 1:length(solution)){
            namenod <- paste('Taxlevel',i,sep='_')
            if (length(taxa) == 0){  # display all if no genera specified
                tmp <- solution[[i]]$result[, "MRCA", drop=FALSE]
                allnodes[[namenod]] <- tmp
            } else {  # display specific genera if specified
                tmp <- solution[[i]]$result[taxa, "MRCA", drop=FALSE]  # display invading tips/species
                allnodes[[namenod]] <- tmp
            }
        }
    } else {
        if (taxlevels > length(solution)){
            stop('Requested taxonomic level not available (less levels specified as analysis input)!')
        }
        #for (j in 1:length(taxlevels)) {
            namenod <- paste('Taxlevel', taxlevels, sep='_')
            if (length(taxa) == 0){  # display all if no genera specified
                tmp <- solution[[taxlevels]]$result[, "MRCA", drop=FALSE]
                allnodes[[namenod]] <- tmp
            } else {  # display specific genera if specified
                tmp <- solution[[taxlevels]]$result[taxa, "MRCA", drop=FALSE]  # display invading tips/species
                allnodes[[namenod]] <- tmp
            }
        #}
    }
    allnodes
}
