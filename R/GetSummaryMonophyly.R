# get summary only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetSummaryMonophyly <-
function(solution, taxlevels='ALL') {
    Allsummaries <- list()
    if (taxlevels!='ALL' & class(taxlevels)!='numeric') {
	stop("taxlevels must be either 'ALL' or numeric!")
    }
    if (taxlevels=='ALL') {
        for (i in 1:length(solution)){
            namesum <- paste('Taxlevel',i,sep='_')
            tmp <- (solution[[i]]$summary)
            Allsummaries[[namesum]] <- tmp
        }
    } else {
        if (taxlevels > length(solution)){
            stop('Requested taxonomic level not available (less levels specified as analysis input)!')
        }
        #for (j in 1:length(taxlevels)) {
            namesum <- paste('Taxlevel', taxlevels, sep='_')
            tmp <- (solution[[taxlevels]]$summary)
            Allsummaries[[namesum]] <- tmp
        #}
    }
    Allsummaries
}
