#get list of intruders (genera) only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetIntruderTaxa <-
function(solution, taxa=NULL, taxlevels='ALL') {
    alltaxa <- list()
    if (taxlevels!='ALL' & class(taxlevels)!='numeric') {
	stop("taxlevels must be either 'ALL' or numeric!")
    }
    if (taxlevels=='ALL') {
        for (i in 1:length(solution)){
            nametax <- paste('Taxlevel',i,sep='_')
            if (length(taxa) == 0){  # display all if no taxon specified
                tmp <- solution[[i]]$IntruderTaxa
                alltaxa[[nametax]] <- tmp
            } else {  # display specific taxon if specified
                alltaxa2 <- list()
                for (i in 1:length(taxa)) {  # loop to go through vector of taxon names
                    nametax2 <- taxa[i]  # display name of invaded taxon first
                    tmp <- solution[[i]]$IntruderTaxa[[taxa[i]]]  # display invading taxa
                    alltaxa2[[nametax2]] <- tmp
                }
                alltaxa[[nametax]] <- alltaxa2
            }
        }
    } else {
        if (taxlevels > length(solution)){
            stop('Requested taxonomic level not available (less levels specified as analysis input)!')
        }
                #for (j in 1:length(taxlevels)) {
            nametax <- paste('Taxlevel', taxlevels, sep='_')
            if (length(taxa) == 0){  # display all if no taxa specified
                tmp <- solution[[taxlevels]]$IntruderTaxa
                alltaxa[[nametax]] <- tmp
            } else {  # display specific taxa if specified
                alltaxa2 <- list()
                for (i in 1:length(taxa)) {  # loop to go through vector of taxon names
                    nametax2 <- taxa[i]  # display name of invaded taxon first
                    tmp <- solution[[taxlevels]]$IntruderTaxa[[taxa[i]]]  # display invading taxa
                    alltaxa2[[nametax2]] <- tmp
                }
                alltaxa[[nametax]] <- alltaxa2
            }
        #}
    }
    alltaxa
}
