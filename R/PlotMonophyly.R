# plot functions to visualize the outputs of AssessMonophyly
# written by Orlando Schwery 2015

PlotMonophyly <-
function(solution, tree, type='momophyly', ladderize=TRUE, PDF=FALSE, PDF_filename='Monophylyplot.pdf', mono.colour='gwp', tax.colour='rainbow', intrud.colour='rainbow', ...) {
    
    if (ladderize==TRUE) {  # ladderizes the tree before starting, if specified
        tree <- ladderize(tree)
    }

# Reconstruct Monophyly
    
    if (type=='monophyly' | type=='monoVStax' | type=='intruders') {  # for the plot variants that require monophyly status reconstructed
        mono.tree <- tree  # assigns tree specifically to reconstruction to avoid conflicts
        
        # assign numbers to tip status
        tipdata <- as.character(solution$TipStates[, "Status"])  # extracting monophyly status of tips from solution
        tipdata[tipdata == "Monophyletic"] <- 1  # number-coding monophyly status
        tipdata[tipdata == "Non-Monophyletic"] <- 2  # number-coding monophyly status
        tipdata[tipdata == "Intruder"] <- 3  # number-coding monophyly status
        tipdata <- as.numeric(tipdata)
         # run reco model
        monophyly.reco <- fastAnc(mono.tree, tipdata, vars=FALSE, CI=FALSE)
    
        edgestates <- c()  # empty vector for reconstructed states at edges
        for (i in 1:length(mono.tree$edge[, 2])) {
            if (mono.tree$edge[i, 2] > length(mono.tree$tip.label)){  # if internal edge
                edgestates.i <- as.vector(monophyly.reco[mono.tree$edge[i, 2] - (length(mono.tree$tip.label))])  # use reconstructed edge state
            }else{  # if terminal branch
                edgestates.i <- as.vector(tipdata[mono.tree$edge[i,2]])  # use tip state
            }
            edgestates <- c(edgestates, edgestates.i)  # add to vector
        }
        
        mono.tree$edge <- cbind(mono.tree$edge, round(as.numeric(edgestates), digits=0))  # round values to full digits (i.e. initial states) and add to edges of tree
    }
    
##########################
    # reconstruct intruders in monophyly
    if (type=='intruders') {
        int.tree <- tree  # assigns tree specifically to reconstruction to avoid conflicts
        
        # assign numbers to tip status
        tipdataI <- as.character(solution$TipStates[, "Status"])  # extracting monophyly status of tips from solution
        tipdataI[tipdataI == "Monophyletic"] <- 1  # number-coding monophyly status
        tipdataI[tipdataI == "Non-Monophyletic"] <- 2  # number-coding monophyly status
        tipdataI[tipdataI == "Intruder"] <- 3  # number-coding monophyly status
        tipdataI <- as.factor(tipdataI)
        
        #assign numbers to tip genus
        tipdataII <- as.character(solution$TipStates[, "Genus"])  #vector with genus names
        taxai <- c()
        for (i in 1:length(tipdataI)){
            if (tipdataI[i] == 3){  # list taxa of tips classified as intruders
                taxai <- c(taxai, tipdataII[i])
            }
        }
        taxaI <- as.vector(unique(taxai))  # create vector of intruder taxa (without doubles)
        
        taxaII <- c()
        for (i in 1:length(taxaI)){
            taxaII[i] <- i-1  # create vector representing numbers to intruding taxon
        }
        taxaIII <- matrix(c(taxaI, taxaII), nrow=length(taxaI)) #create matrix consisting of intruder taxon names and their number
        
        tipdataIII <- c()  # create empty vector for tip states to be used in reconstruction
        for (i in 1:length(tipdataI)) {
            if (tipdataI[i] == 3){
                tipdataIII <- c(tipdataIII,(as.numeric(tipdataI[i]) + as.numeric(taxaIII[(which(taxaIII[, 1] == tipdataII[i])), 2])))  # if intruder, add previously assigned taxon number (leading intruders to be coded from 3 onwards)
            } else {
                tipdataIII <- c(tipdataIII, as.numeric(tipdataI[i]))  # if not intruder, keep number 1 or 2 respectively
            }
        }
        tipdataIII <- as.numeric(tipdataIII)
        
        # run reco model
        int.reco <- fastAnc(int.tree, tipdataIII, vars=FALSE, CI=FALSE)

        edgestatesI <- c()  # empty vector for reconstructed states at edges
        for (i in 1:length(int.tree$edge[, 2])) {
            if (int.tree$edge[i, 2] > length(int.tree$tip.label)){  # if internal edge
                edgestatesI.i <- as.vector(int.reco[int.tree$edge[i, 2] - (length(int.tree$tip.label))])  # use reconstructed edge state
            }else{  # if terminal branch
                edgestatesI.i <- as.vector(tipdataIII[int.tree$edge[i,2]])   # use tip state
            }
            edgestatesI <- c(edgestatesI, edgestatesI.i)  # add to vector
        }
        
        edgestates.monoround <- c(round(as.numeric(edgestates), digits=0))
        edgestatesII <- edgestates.monoround
        for (i in 1:length(edgestates)) {
        	if (edgestates.monoround[i]==3) {
        		edgestatesII[i] <- edgestatesI[i]
        	}
        }
        int.tree$edge <- cbind(int.tree$edge, as.numeric(edgestatesII))  # add to edges of tree
    }

##########################
    # reconstruct taxonomy
    if (type=='monoVStax' | type=='taxonomy') {  # for the plot variants that require monophyly status reconstructed
        tax.tree <- tree  # assigns tree specifically to reconstruction to avoid conflicts
        tipdataT <- as.character(solution$TipStates[, "Genus"])  # extract taxon of tip from solution
        taxaT <- as.vector(unique(solution$TipStates[, "Genus"])) # vector with taxon names (no doubles)
        for (i in 1:length(taxaT)){
            tipdataT[tipdataT == taxaT[i]] <- i  # translate taxon associated to tip into taxon specific number
        }
        tipdataT <- as.numeric(tipdataT)
        
        # run reco model
        taxa.reco <- fastAnc(tax.tree, tipdataT, vars=FALSE, CI=FALSE)

        edgestatesT <- c()  # empty vector for reconstructed states at edges
        for (i in 1:length(tax.tree$edge[, 2])) {
            if (tax.tree$edge[i, 2] > length(tax.tree$tip.label)){  # if internal edge
                edgestatesT.i <- as.vector(taxa.reco[tax.tree$edge[i, 2] - (length(tax.tree$tip.label))])  # use reconstructed edge state
            }else{  # if terminal branch
                edgestatesT.i <- as.vector(tipdataT[tax.tree$edge[i,2]])    # use tip state
            }
            edgestatesT <- c(edgestatesT, edgestatesT.i)  # add to vector
        }
        tax.tree$edge <- cbind(tax.tree$edge, as.numeric(edgestatesT))  # add to edges of tree
    }

################        
    # plotting itself
    # monophyly plot
    if (type=='monophyly' ) {
        if (mono.colour=='gwp') {
            co <- rev(brewer.pal(3, 'PRGn'))  # use predefined colour brewer colours
        } else {
            co <- mono.colour  # use custom colours
        }

        if (PDF==TRUE) {
            pdf(PDF_filename, width=7, height=(length(mono.tree$tip.label)/10))  # create PDF with lenght adjusted to tree size
            plot(mono.tree, edge.col=co[as.numeric(mono.tree$edge[, 3])], cex=0.2, label.offset=3, edge.width=2, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = co[tipdata], cex = 1, adj = 1)  # add tip labels with tip state colour
            dev.off()
        } else {
            plot(mono.tree, edge.col=co[as.numeric(mono.tree$edge[, 3])], cex=0.2, label.offset=3, edge.width=2, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = co[tipdata], cex = 1, adj = 1)  # add tip labels with tip state colour
        }
    }
    
    # taxonomy plot
    if (type=='taxonomy' ) {
        if (tax.colour=='rainbow') {
            coTax <- rainbow(length(taxaT))  # use default taxonomy rainbow colours based on number of taxa
        } else {
            coTax <- tax.colour  # use custom taxonomy colours
        }        
        names(coTax) <- 1:length(taxaT) 

        if (PDF==TRUE) {
            pdf(PDF_filename, width=7, height=(length(tax.tree$tip.label)/10))   # create PDF with lenght adjusted to tree size
            plot(tax.tree, edge.col=coTax[tax.tree$edge[, 3]], cex=0.2, label.offset=3, edge.width=2, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = coTax[tipdataT], cex = 1, adj = 1)  # add tip labels with tip state colour
            dev.off()
        } else {
            plot(tax.tree, edge.col=coTax[tax.tree$edge[, 3]], cex=0.2, label.offset=3, edge.width=2, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = coTax[tipdataT], cex = 1, adj = 1)  # add tip labels with tip state colour
        }
    }

    # intruders plot
    if (type=='intruders') {
         if (intrud.colour=='rainbow') {
            coInt <- c('gray', 'black', rainbow(length(taxaI)))  # use default colours (gray for monophyletic, black for invaded and rainbow colours for invading taxa)
        } else {
            coInt <- intrud.colour  # use custom intruder colours
        }
        names(coInt) <- 1:(length(taxaI) + 2)
   
        if (PDF==TRUE) {
            pdf(PDF_filename, width=7, height=(length(int.tree$tip.label)/10))   # create PDF with lenght adjusted to tree size
            plot(int.tree, edge.col=coInt[int.tree$edge[, 3]], show.tip.label = TRUE, cex=0.2, label.offset=3, edge.width=2, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = coInt[tipdataIII], cex = 1, adj = 1)  # add tip labels with tip state colour
            dev.off()
        } else {     
            plot(int.tree, edge.col=coInt[int.tree$edge[, 3]], show.tip.label = TRUE, cex=0.2, label.offset=3, edge.width=2, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = coInt[tipdataIII], cex = 1, adj = 1)  # add tip labels with tip state colour
        } 
    }
    

       # monophyly vs taxonomy mirror tree plot
    if (type=='monoVStax') {
        if (mono.colour=='gwp') {
            co <- rev(brewer.pal(3, 'PRGn'))  # use default colour brewer colours
        } else {
            co <- mono.colour  # use custom monophyly colour
        }
        if (tax.colour=='rainbow') {
            coTax <- rainbow(length(taxaT))  # use default rainbow colours based on number of taxa
        } else {
            coTax <- tax.colour  # use custom taxa colour
        }        
        names(coTax) <- 1:length(taxaT) 

        if (PDF==TRUE) {
            pdf(PDF_filename, width=14, height=(length(tax.tree$tip.label)/10))   # create PDF with lenght adjusted to tree size
            par(oma=c(1, 1, 1, 1), mar=c(0,0,0,0))  # set up plotting margins
            layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(2,1))  # set up for plotting two trees next to each other
            plot(mono.tree, edge.col=co[as.numeric(mono.tree$edge[, 3])], cex=0.2, adj=0.5, label.offset=15, edge.width=2, no.margin=TRUE, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = co[tipdata], cex = 1, adj = 1)  # add tip labels with tip state colour
            plot(tax.tree, edge.col=coTax[tax.tree$edge[, 3]], show.tip.label = FALSE, cex=0.2, label.offset=(-3), edge.width=2, direction = "leftwards", no.margin=TRUE, ...)  # plot mirrored tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = coTax[tipdataT], cex = 1, adj = 1)  # add tip labels with tip state colour
            dev.off()
        } else {
            par(oma=c(1, 1, 1, 1), mar=c(0,0,0,0))  # set up plotting margins
            layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(2,1))  # set up for plotting two trees next to each other
            plot(mono.tree, edge.col=co[as.numeric(mono.tree$edge[, 3])], cex=0.2, adj=0.5, label.offset=15, edge.width=2, no.margin=TRUE, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = co[tipdata], cex = 1, adj = 1)  # add tip labels with tip state colour
            plot(tax.tree, edge.col=coTax[tax.tree$edge[, 3]], show.tip.label = FALSE, cex=0.2, label.offset=(-3), edge.width=2, direction = "leftwards", no.margin=TRUE, ...)  # plot mirrored tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = coTax[tipdataT], cex = 1, adj = 1)  # add tip labels with tip state colour
        }

    }
    

}
