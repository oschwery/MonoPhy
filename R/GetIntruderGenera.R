#get list of intruders (genera) only displayed out of the outputlist of AssessMonophyly
# written by Orlando Schwery 2015

GetIntruderGenera <-
function(solution, genera=NULL) {
    if (length(genera) == 0){  # display all if no genera specified
    print(solution$Genera)
    } else {  # display specific genera if specified
        for (i in 1:length(genera)) {  # loop to go through vector of genus names
            writeLines(genera[i])  # display name of invaded genus first
            print(solution$Genera[[genera[i]]])  # display invading genera
            writeLines("\n")  # add empty line for better readability
        }
    }
}
