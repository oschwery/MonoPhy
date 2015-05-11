#get list of intruders (tips) only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetIntruderTips <-
function(solution, taxa=NULL, taxlevels='ALL') {
    alltips <- list()
    if (taxlevels=='ALL') {
        for (i in 1:length(solution)){
            nametip <- paste('Taxlevel',i,sep='_')
            if (length(taxa) == 0){  # display all if no taxon specified
                tmp <- solution[[i]]$Species
                alltips[[nametip]] <- tmp
            } else {  # display specific taxon if specified
                alltips2 <- list()
                for (i in 1:length(taxa)) {  # loop to go through vector of taxon names
                    nametip2 <- taxa[i]  # display name of invaded taxon first
                    tmp <- solution[[i]]$Species[[taxa[i]]]  # display invading taxa
                    alltips2[[nametip2]] <- tmp
                }
                alltips[[nametip]] <- alltips2
            }
        }
    } else {
        if (taxlevels > length(solution)){
            stop('Requested taxonomic level not available (less levels specified as analysis input)!')
        }
                #for (j in 1:length(taxlevels)) {
            nametip <- paste('Taxlevel', taxlevels, sep='_')
            if (length(taxa) == 0){  # display all if no genera specified
                tmp <- solution[[taxlevels]]$Species
                alltips[[nametip]] <- tmp
            } else {  # display specific genera if specified
                alltips2 <- list()
                for (i in 1:length(taxa)) {  # loop to go through vector of taxon names
                    nametip2 <- taxa[i]  # display name of invaded taxon first
                    tmp <- solution[[taxlevels]]$Species[[taxa[i]]]  # display invading taxa
                    alltips2[[nametip2]] <- tmp
                }
                alltips[[nametip]] <- alltips2
            }
        #}
    }
    alltips
}
