# plot functions to visualize the outputs of AssessMonophyly
# written by Orlando Schwery 2015

PlotMonophyly <-
function(solution, tree, taxlevels=1, plot.type='monophyly', monocoll=FALSE, ladderize=TRUE, PDF=FALSE, PDF_filename='Monophylyplot.pdf', PDF_width='auto', PDF_height='auto', mono.colour='PRGn', tax.colour='rainbow', intrud.colour='rainbow', edge.width=3, cex=0.2, type='phylogram', ...) {
    
    if (taxlevels == 'ALL' | class(taxlevels)!='numeric') {
		stop("taxlevels must be numeric (also 'ALL' is not an option for plotting)!")
    }    
    if (taxlevels > length(solution)) {
		stop('Requested taxonomic level not available (less levels specified as analysis input)!')
    }
    if (plot.type!='monophyly' & plot.type!='monoVStax' & plot.type!='intruders' & plot.type!='taxonomy') {
		stop('Invalid plot.type!')
    }
    if (type == "radial") {
		stop("Type 'radial' is currently not supported!")
	}
	if (ladderize==TRUE) {  # ladderizes the tree before starting, if specified
        tree <- ladderize(tree)
    }
    tip.states <- solution[[taxlevels]]$TipStates
    tip.states <- tip.states[match(tree$tip.label, tip.states$Tip),]
	row.names(tip.states) <- 1:nrow(tip.states)
	
    if (monocoll== TRUE) {
		tip.states$Tip <- as.character(tip.states$Tip)
		tip.states$Taxon <- as.character(tip.states$Taxon)
		tip.states$Status <- as.character(tip.states$Status)		

		alltaxa <- as.vector(unique(tip.states[, "Taxon"]))
		keeptips <- c()
		colltipse <- c()
		for (icoll in 1:length(alltaxa)) {
		    if (solution[[taxlevels]]$result[alltaxa[icoll], "Monophyly"] == "Yes") {
			matchtips <- which(tip.states[, "Taxon"] == alltaxa[icoll])
			colltips <- c(matchtips[2:length(matchtips)])
			colltipse <- c(colltipse,colltips)
			keeptip <- c(matchtips[1])
			keeptips <- c(keeptips, keeptip)
		    }
		}
		
		tip.states.temp <- tip.states
		
		for (icollstates in 1:length(colltipse)) {
		    tip.states.temp <- tip.states.temp[!(tip.states.temp$Tip == tree$tip.label[colltipse[icollstates]]),]
		}
		row.names(tip.states.temp) <- 1:nrow(tip.states.temp)

		for (inameadj in 1:length(tip.states.temp$Tip)) {
			if (solution[[taxlevels]]$result[tip.states.temp$Taxon[inameadj], "Monophyly"] == "Yes") {
				tip.states.temp$Tip[inameadj] <- tip.states.temp$Taxon[inameadj]
			}
		}
		
		newlabels <- tree$tip.label
		for (itips in 1:length(keeptips)) {
		    newlabels[keeptips[itips]] <- as.character(tip.states[keeptips[itips],"Taxon"])
		}
		tree$tip.label <- newlabels
		tree <- drop.tip(tree,colltipse)
		tip.states <- tip.states.temp
    }
    
# Reconstruct Monophyly
    
    if (plot.type=='monophyly' | plot.type=='monoVStax' | plot.type=='intruders') {  # for the plot variants that require monophyly status reconstructed
        mono.tree <- tree  # assigns tree specifically to reconstruction to avoid conflicts
        
        # assign numbers to tip status
        tipdata <- as.character(tip.states[, "Status"])  # extracting monophyly status of tips from solution
        tipdata[tipdata == "Monophyletic"] <- 2  # number-coding monophyly status
        tipdata[tipdata == "Non-Monophyletic"] <- 3  # number-coding monophyly status
        tipdata[tipdata == "Intruder"] <- 4  # number-coding monophyly status
        tipdata[tipdata == "Outlier"] <- 4  # number-coding monophyly status
		tipdata[tipdata == "unknown"] <- 1  # number-coding monophyly status
		tipdata <- as.numeric(tipdata)
		names(tipdata) <- tip.states[, "Tip"]
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
    if (plot.type=='intruders') {
        int.tree <- tree  # assigns tree specifically to reconstruction to avoid conflicts
        
        # assign numbers to tip status
        tipdataI <- as.character(tip.states[, "Status"])  # extracting monophyly status of tips from solution
        tipdataI[tipdataI == "Monophyletic"] <- 2  # number-coding monophyly status
        tipdataI[tipdataI == "Non-Monophyletic"] <- 3  # number-coding monophyly status
        tipdataI[tipdataI == "Intruder"] <- 4  # number-coding monophyly status
        tipdataI[tipdataI == "Outlier"] <- 4  # number-coding monophyly status
	tipdataI[tipdataI == "unknown"] <- 1  # number-coding monophyly status
        tipdataI <- as.numeric(tipdataI)
        
        #assign numbers to tip genus
        tipdataII <- as.character(tip.states[, "Taxon"])  #vector with genus names
        taxai <- c()
        for (i in 1:length(tipdataI)){
            if (tipdataI[i] == 4){  # list taxa of tips classified as intruders
                taxai <- c(taxai, tipdataII[i])
            }
        }
        if (length(taxai)==0) {
            stop('No intruders present to be plotted! Get on with it!')
        }
	taxaI <- as.vector(unique(taxai))  # create vector of intruder taxa (without doubles)
        
        taxaII <- c()
        for (i in 1:length(taxaI)){
            taxaII[i] <- i-1  # create vector representing numbers to intruding taxon
        }
        taxaIII <- matrix(c(taxaI, taxaII), nrow=length(taxaI)) #create matrix consisting of intruder taxon names and their number
        
        tipdataIII <- c()  # create empty vector for tip states to be used in reconstruction
        for (i in 1:length(tipdataI)) {
            if (tipdataI[i] == 4){
                tipdataIII <- c(tipdataIII,(as.numeric(tipdataI[i]) + as.numeric(taxaIII[(which(taxaIII[, 1] == tipdataII[i])), 2])))  # if intruder, add previously assigned taxon number (leading intruders to be coded from 4 onwards)
            } else {
                tipdataIII <- c(tipdataIII, as.numeric(tipdataI[i]))  # if not intruder, keep number 1, 2 or 3 respectively
            }
        }
	tipdataIII <- as.numeric(tipdataIII)
	names(tipdataIII) <- tip.states[, "Tip"]
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
        	if (edgestates.monoround[i]==4) {
        		edgestatesII[i] <- edgestatesI[i]
        	}
        }
	int.tree$edge <- cbind(int.tree$edge, as.numeric(edgestatesII))  # add to edges of tree
    }

##########################
    # reconstruct taxonomy
    if (plot.type=='monoVStax' | plot.type=='taxonomy') {  # for the plot variants that require monophyly status reconstructed
        tax.tree <- tree  # assigns tree specifically to reconstruction to avoid conflicts
        tipdataT <- as.character(tip.states[, "Taxon"])  # extract taxon of tip from solution
        taxaT <- as.vector(unique(tip.states[, "Taxon"])) # vector with taxon names (no doubles)
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
    if (PDF==TRUE) {
		if (PDF_width=='auto') {
			if (type =="fan") {
				pdf_width <- 2.5*sqrt(length(tree$tip.label))/pi  # create PDF with width adjusted to tree size (square shaped for round trees)
			} else {
				pdf_width <- (9-(3-(3*(2^-(length(tree$tip.label)/100)))))  # create PDF with width adjusted to tree size (rectangular for straight trees)
			}
		} else {
			pdf_width <- PDF_width
		}
		if (PDF_height=='auto') {
			if (type =="fan") {
				pdf_height <- 2.5*sqrt(length(tree$tip.label))/pi  # create PDF with lenght adjusted to tree size (square shaped for round trees)
			} else {
				pdf_height <- (length(tree$tip.label)/10)  # create PDF with lenght adjusted to tree size (rectangular for straight trees)
			}
		} else {
			pdf_height <- PDF_height
		}
		if (pdf_height > 200) {
			pdf_height <- 200
			print('Warning: pdf_height too large, capped to 200in. If output not satisfying, consider different plotting type or tree slicing.')
		}
		if (pdf_width > 200) {
			pdf_width <- 200
			print('Warning: pdf_width too large, capped to 200in. If output not satisfying, consider different plotting type or tree slicing.')
		}
	}
    # monophyly plot
    if (plot.type=='monophyly' ) {
        if (mono.colour=='PRGn') {  # use colours from colourblind friendly palettes of RColorBrewer
	    co <- c('gray', '#5aae61', '#c2a5cf', '#762a83') #green and purple 
	} else if (mono.colour=='RdBu') {
	    co <- c('gray', '#4393c3', '#f4a582', '#b2182b') #red and blue
	} else if (mono.colour=='PuOr') {
	    co <- c('gray', '#542788', '#fdb863', '#b35806') #orange and purple
	} else if (mono.colour=='PiYG') {
	    co <- c('gray', '#4d9221','#f1b6da', '#8e0152') #green and pink
	} else if (mono.colour=='BrBG') {
	    co <- c('gray', '#35978f','#dfc27d', '#543005') #petrol and brown
	} else {
            co <- mono.colour  # use custom colours
        }

        if (PDF==TRUE) {
	    pdf(PDF_filename, width=pdf_width, height=pdf_height)  # create PDF frame
            plot(mono.tree, edge.col=co[as.numeric(mono.tree$edge[, 3])], cex=cex, label.offset=3-(2-(2.5*(2^-(length(mono.tree$tip.label)/100)))), edge.width=edge.width, type=type, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = co[tipdata], cex = 1, adj = 1)  # add tip labels with tip state colour
            dev.off()
        } else {
            plot(mono.tree, edge.col=co[as.numeric(mono.tree$edge[, 3])], cex=cex, label.offset=3, edge.width=edge.width, type=type, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = co[tipdata], cex = 1, adj = 1)  # add tip labels with tip state colour
        }
    }
    
    # taxonomy plot
    if (plot.type=='taxonomy' ) {
        if (tax.colour=='rainbow') {
            coTax <- rainbow(length(taxaT))  # use default taxonomy rainbow colours based on number of taxa
        } else {
            coTax <- tax.colour  # use custom taxonomy colours
        }        
        names(coTax) <- 1:length(taxaT) 

        if (PDF==TRUE) {
	    pdf(PDF_filename, width=pdf_width, height=pdf_height)  # create PDF frame
            plot(tax.tree, edge.col=coTax[as.numeric(tax.tree$edge[, 3])], cex=cex, label.offset=3-(2-(2.5*(2^-(length(tax.tree$tip.label)/100)))), edge.width=edge.width, type=type, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = coTax[tipdataT], cex = 1, adj = 1)  # add tip labels with tip state colour
            dev.off()
        } else {
            plot(tax.tree, edge.col=coTax[as.numeric(tax.tree$edge[, 3])], cex=cex, label.offset=3, edge.width=edge.width, type=type, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = coTax[tipdataT], cex = 1, adj = 1)  # add tip labels with tip state colour
        }
    }

    # intruders plot
    if (plot.type=='intruders') {
	    if (intrud.colour=='rainbow') {
	        coInt <- c('gray82', 'gray50', 'black', rainbow(length(taxaI)))  # use default colours (gray for monophyletic, black for invaded and rainbow colours for invading taxa)
	    } else {
	        coInt <- intrud.colour  # use custom intruder colours
	    }
        names(coInt) <- 1:(length(taxaI) + 2)
   
        if (PDF==TRUE) {
	    pdf(PDF_filename, width=pdf_width, height=pdf_height)  # create PDF frame
	    plot(int.tree, edge.col=coInt[as.numeric(int.tree$edge[, 3])], show.tip.label = TRUE, cex=cex, label.offset=3-(2-(2.5*(2^-(length(int.tree$tip.label)/100)))), edge.width=edge.width, type=type, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = coInt[tipdataIII], cex = 1, adj = 1)  # add tip labels with tip state colour
            dev.off()
        } else {     
            plot(int.tree, edge.col=coInt[as.numeric(int.tree$edge[, 3])], show.tip.label = TRUE, cex=cex, label.offset=3, edge.width=edge.width, type=type, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = coInt[tipdataIII], cex = 1, adj = 1)  # add tip labels with tip state colour
        } 
    }
    

       # monophyly vs taxonomy mirror tree plot
    if (plot.type=='monoVStax') {
        if (type == 'fan' | type == 'unrooted') {
	    stop("Phylogeny types 'fan' and 'unrooted' aren't suited for plot.type 'monoVStax', use types 'phylogram' or 'cladogram' instead!")
        }
        if (mono.colour=='PRGn') {  # use colours from colourblind friendly palettes of RColorBrewer
	    co <- c('gray', '#5aae61','#c2a5cf', '#762a83') #green and purple 
	} else if (mono.colour=='RdBu') {
	    co <- c('gray', '#4393c3','#f4a582', '#b2182b') #red and blue
	} else if (mono.colour=='PuOr') {
	    co <- c('gray', '#542788','#fdb863', '#b35806') #orange and purple
	} else if (mono.colour=='PiYG') {
	    co <- c('gray', '#4d9221','#f1b6da', '#8e0152') #green and pink
	} else if (mono.colour=='BrBG') {
	    co <- c('gray', '#35978f','#dfc27d', '#543005') #petrol and brown
	} else {
            co <- mono.colour  # use custom colours
        }
        if (tax.colour=='rainbow') {
            coTax <- rainbow(length(taxaT))  # use default rainbow colours based on number of taxa
        } else {
            coTax <- tax.colour  # use custom taxa colour
        }        
        names(coTax) <- 1:length(taxaT) 

        if (PDF==TRUE) {
            pdf(PDF_filename, width=pdf_width, height=pdf_height)   # create PDF with lenght adjusted to tree size
            par(oma=c(1, 1, 1, 1), mar=c(0,0,0,0))  # set up plotting margins
            layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(1.5,1))  # set up for plotting two trees next to each other
            plot(mono.tree, edge.col=co[as.numeric(mono.tree$edge[, 3])], cex=cex, adj=0.5, label.offset=200000/(length(mono.tree$tip.label)^2), edge.width=edge.width, no.margin=TRUE, type=type, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = co[tipdata], cex = 1, adj = 0.5)  # add tip labels with tip state colour
            plot(tax.tree, edge.col=coTax[as.numeric(tax.tree$edge[, 3])], show.tip.label = FALSE, cex=cex, label.offset=0, edge.width=edge.width, direction = "leftwards", no.margin=TRUE, type=type, ...)  # plot mirrored tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = coTax[tipdataT], cex = 1, adj = 0.5)  # add tip labels with tip state colour
            dev.off()
        } else {
            par(oma=c(1, 1, 1, 1), mar=c(0,0,0,0))  # set up plotting margins
            layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(1.5,1))  # set up for plotting two trees next to each other
            plot(mono.tree, edge.col=co[as.numeric(mono.tree$edge[, 3])], cex=cex, adj=0.5, label.offset=200000/(length(mono.tree$tip.label)^2), edge.width=edge.width, no.margin=TRUE, type=type, ...)  # plot tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = co[tipdata], cex = 1, adj = 0.5)  # add tip labels with tip state colour
            plot(tax.tree, edge.col=coTax[as.numeric(tax.tree$edge[, 3])], show.tip.label = FALSE, cex=cex, label.offset=0, edge.width=edge.width, direction = "leftwards", no.margin=TRUE, type=type, ...)  # plot mirrored tree with edge colours according to reconstruction
            tiplabels(pch = 22, bg = coTax[tipdataT], cex = 1, adj = 0.5)  # add tip labels with tip state colour
        }
    }
}