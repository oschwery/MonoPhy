#get list of outliers (genera) only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetOutlierTaxa <-
function(solution, taxlevels='ALL') {
    alltaxa <- list()
    if (taxlevels!='ALL' & class(taxlevels)!='numeric') {
	stop("taxlevels must be either 'ALL' or numeric!")
    }
    if (taxlevels=='ALL') {
        for (i in 1:length(solution)){
            nametax <- paste('Taxlevel',i,sep='_')
            tmp <- solution[[i]]$OutlierTaxa
            alltaxa[[nametax]] <- tmp
        }
    } else {
        if (taxlevels > length(solution)){
            stop('Requested taxonomic level not available (less levels specified as analysis input)!')
        }
                #for (j in 1:length(taxlevels)) {
            nametax <- paste('Taxlevel', taxlevels, sep='_')
            tmp <- solution[[taxlevels]]$OutlierTaxa
            alltaxa[[nametax]] <- tmp
        #}
    }
    alltaxa
}
