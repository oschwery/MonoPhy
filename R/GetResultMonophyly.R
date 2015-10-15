# get result table only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetResultMonophyly <-
function(solution, taxlevels='ALL') {
    Allresults <- list()  # create empty list to be filled
    if (taxlevels != 'ALL' & class(taxlevels) != 'numeric') {  # test format of taxlevels argument and display error if format is wrong
	stop("taxlevels must be either 'ALL' or numeric!")
    }
    if (taxlevels == 'ALL') {  # if all taxlevels are looked for
        for (i in 1:length(solution)) {  # loop through all taxlevels
            nameres <- paste('Taxlevel', i, sep='_')  # create namelabel for current taxlevel
            tmp <- (solution[[i]]$result)  # extract result object from solution
            Allresults[[nameres]] <- tmp  # add extracted table as sub-object to output list and label it with the appropriate taxlevel nr.
        }
    } else {  # if only a specific taxlevel is requested
        if (taxlevels > length(solution)) {  # test whether requested taxlevel is among available ones and display error if not
            stop('Requested taxonomic level not available (less levels specified as analysis input)!')
        }
        nameres <- paste('Taxlevel', taxlevels, sep='_')  # create namelabel for current taxlevel
        tmp <- (solution[[taxlevels]]$result)  # extract result object from solution
        Allresults[[nameres]] <- tmp  # add extracted list as sub-object to output list and label it with the appropriate taxlevel nr.
    }
    Allresults  # export output list
}
