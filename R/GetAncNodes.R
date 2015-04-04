#get list of MRCA nodes of desired genera only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetAncNodes <-
function(solution, genera=NULL) {
    if (length(genera) == 0){  # display all if no genera specified
    print(solution$result[, "MRCA", drop=FALSE])
    } else {  # display specific genera if specified
        print(solution$result[genera, "MRCA", drop=FALSE])  # display invading tips/species
    }
}
