#get list of intruders (tips) only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetIntruderTips <-
function(solution, genera=NULL) {
    if (length(genera) == 0){  # display all if no genera specified
    print(solution$Species)
    } else {  # display specific genera if specified
        for (i in 1:length(genera)) {  # loop to go through vector of genus names
            writeLines(genera[i])  # display name of invaded genus first
            print(solution$Species[[genera[i]]])  # display invading tips/species
            writeLines("\n")  # add empty line for better readability
        }
    }
}
