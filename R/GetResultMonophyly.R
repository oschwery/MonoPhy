# get result table only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetResultMonophyly <-
function(solution, taxlevels='ALL') {
    Allresults <- list()
    if (taxlevels=='ALL') {
        for (i in 1:length(solution)){
            nameres <- paste('Taxlevel',i,sep='_')
            tmp <- (solution[[i]]$result)
            Allresults[[nameres]] <- tmp
        }
    } else {
        if (taxlevels > length(solution)){
            stop('Requested taxonomic level not available (less levels specified as analysis input)!')
        }
        #for (j in 1:length(taxlevels)) {
            nameres <- paste('Taxlevel', taxlevels, sep='_')
            tmp <- (solution[[taxlevels]]$result)
            Allresults[[nameres]] <- tmp
        #}
    }
    Allresults
}
