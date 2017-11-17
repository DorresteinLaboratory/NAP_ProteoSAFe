#!/PATHCONDAENV/envs/nap/bin/Rscript

library(hwriter)
library(fingerprint)
library(rcdk)
library(ggplot2)
library(dynamicTreeCut)
library(metfRag)
library(rjson)
suppressPackageStartupMessages(library(optparse))


  load("split_data/tabgnps.rda")
  load("lid_res/lid.rda")
  load("merge_fusion/fusion.rda")
  load("merge_consensus/consensus.rda")

  jtab <- fromJSON(file='split_data/config.json')
  args <- commandArgs(TRUE)

  PATH <- "/PATHTOOLSFOLDER/nap_ccms2/Snap/"
  opt <- list()
  if(PATH!=""){
	opt$path <- PATH
  }  

  source(paste0(PATH, "code/load.db.R")) 

#  if(jtab$database=="PubChem" | jtab$database=="PubChem[Only]") {
#	  TYPE <- "PubChem"
#  } else {
	  TYPE <- "local"
	  if(jtab$database=="none" & jtab$udb!="0") {
		database = fread("split_data/udb.txt")
	  } else if(jtab$database!="none" & jtab$udb!="0") {
			dbs <- strsplit(jtab$database, ",")[[1]]
			if(any(grepl("PubChem", dbs))) {
				dbs <- dbs[-grep("PubChem", dbs)]
			}

			database2 = fread("split_data/udb.txt")

			if(length(dbs)) {
				database <- lapply(dbs, function(x) load.db("/PATHTOOLSFOLDER/nap_ccms2/", x))  
				database <- do.call(rbind, database) 
				database = rbind(database, database2)
			}  else{
				database = database2
			}

	  } else {
			dbs <- strsplit(jtab$database, ",")[[1]]
			database <- lapply(dbs, function(x) load.db("/PATHTOOLSFOLDER/nap_ccms2/", x))  
			database <- do.call(rbind, database) 
	  }
#
#   }
   database <- as.data.frame(database) 
   database[, 1] <- as.numeric(database[,1]) 
   if(sum(duplicated(database[,4]))) database <- database[-which(duplicated(database[,4])),]

   dtab <- tabgnps1[,c("cluster.index", "parent.mass", "RTMean", "LibraryID")]

   dtab <- matrix(apply(dtab, 2, unlist), ncol=4)
   dtab  <- gsub("\"", "", dtab)
   dtab <- cbind(dtab, matrix("", nrow=nrow(dtab), ncol=6))
   colnames(dtab) <- c("cluster.index", "parent.mass", "RTMean", "LibraryID", "MetFragID", "MetFragSC", "FusionID", "FusionSC", "ConsensusID", "ConsensusSC" )
   dtab <- gsub("\\s", "", dtab) 

   mlist <- list()
   if (length(jtab[['PLOT']])) {
	dir.create("metfrag_out")
	dir.create("metfrag_fig_out")
   }

	# Make sure that targL2 has the same length and contains all info	
	for(n in 1:length(targL)) {
		if( n <= length(targL2) ) {
			if(!is.null(targL[[n]]) & is.null(targL2[[n]]) ) {
				targL2[[n]] <- targL[[n]]
			} 
		} else {
		 	targL2 <- c(targL2, targL[n])
		}
	}

	# This block generates one summary table for each
	# metfrag result 
	for(i in 1:nrow(tabgnps1)) {
	  if(is.null(lid[[i]])) next
	  targtemp <- lid[[i]]
	  if(is.null(ncol(targtemp))) {
	    n <- names(targtemp)
	    targtemp <- matrix(targtemp, nrow=1)
	    colnames(targtemp) <- n
	  }
	  
	  if(!sum(grepl("Identifier", colnames(targtemp)))) next 
	  if(targtemp[1,"Score"]=="0.0" | targtemp[1,"Score"]=="NA") next
	  
	  if(!is.null(targL2[[i]])) {
	    if(sum(grepl("fusion", colnames(targL2[[i]])))) {
	      if(sum(colnames(targL2[[i]])=="fusion") & sum(colnames(targL2[[i]])=="fusion2")) {
		m <- matrix(cbind(matrix(targtemp[,c("Identifier", "MonoisotopicMass", "Score", "NoExplPeaks")], ncol=4),
				  matrix(targL2[[i]][, c("fusion", "fusion2")], ncol=2)), ncol=6)
		m <- matrix(apply(m, 2, unlist), ncol=6)
		colnames(m) <- c("Identifier", "MonoisotopicMass", "Score", "NoExplPeaks", "fusion", "fusion2")
	      } else if (sum(colnames(targL2[[i]])=="fusion") & !sum(colnames(targL2[[i]])=="fusion2")) {
		m <- matrix(cbind(matrix(targtemp[,c("Identifier", "MonoisotopicMass", "Score", "NoExplPeaks")], ncol=4),
				  targL2[[i]][,"fusion"]), ncol=5)
		m <- matrix(apply(m, 2, unlist), ncol=5)
		colnames(m) <- c("Identifier", "MonoisotopicMass", "Score", "NoExplPeaks", "fusion")
	      } else if (!sum(colnames(targL2[[i]])=="fusion") & sum(colnames(targL2[[i]])=="fusion2")) {
		m <- matrix(cbind(matrix(targtemp[,c("Identifier", "MonoisotopicMass", "Score", "NoExplPeaks")], ncol=4),
				  targL2[[i]][,"fusion2"]), ncol=5)
		m <- matrix(apply(m, 2, unlist), ncol=5)
		colnames(m) <- c("Identifier", "MonoisotopicMass", "Score", "NoExplPeaks", "fusion2")
	      } 
	    } else {
	      m <- matrix(targtemp[,c("Identifier", "MonoisotopicMass", "Score",  "NoExplPeaks")], ncol=4)
	      m <- matrix(apply(m, 2, unlist), ncol=4)
	      colnames(m) <- c("Identifier", "MonoisotopicMass", "Score", "NoExplPeaks")
	    } 
	  } else {
	    m <- matrix(targtemp[,c("Identifier", "MonoisotopicMass", "Score", "NoExplPeaks")], ncol=4)
	    m <- matrix(apply(m, 2, unlist), ncol=4)
	    colnames(m) <- c("Identifier", "MonoisotopicMass", "Score", "NoExplPeaks")
	  }

	  if(nrow(m)>10) n <- 10 else n <- nrow(m)
	  
	  if(nrow(m)>2) {
		  mols <- parse.smiles(as.character(lid[[i]][,"SMILES"]))
		  fp <- lapply(mols, function(x) try(get.fingerprint(x), TRUE))
		  perr <- unlist(lapply(fp, function(x) any(is(x)!="fingerprint")))
		  if(sum(perr)) {
			fpm <- matrix(0, nrow=length(mols), ncol=length(mols))
		  	fpmtmp <- fp.sim.matrix(fp[-which(perr)])
			fpm[-which(perr),-which(perr)] <- fpmtmp 
		  } else {
		  	fpm <- fp.sim.matrix(fp)
		  }

                  if(nrow(fpm)<11) mds <- cmdscale( as.dist(1-fpm), k=nrow(fpm)-1) else mds <- cmdscale( as.dist(1-fpm), k=10)
		  if(ncol(mds)==0) {
			 mtmp <- as.matrix(database[match(m[,1], database[,4]), c("superclass_name", "class_name")], nrow=nrow(m))
			 m <- cbind(matrix(m[,1], nrow=nrow(m)), matrix('', ncol=1, nrow=nrow(m)), mtmp,  matrix(m[,-1], nrow=nrow(m)) )
			 if(sum(grepl("Score", colnames(m)))) {
				sup <- unique(m[, "superclass_name"]) 
				sup <- paste(sup, collapse="#") 
				id <- unique(m[, "Identifier"])
				id <- paste(id, collapse="#") 
				dtab[i,"MetFragID"] <- id
				dtab[i,"MetFragSC"] <- sup 
			 }
			 if(any(colnames(m)=="fusion" )) {
				sup <- unique(m[, "superclass_name"]) 
				sup <- paste(sup, collapse="#") 
				id <- unique(m[, "Identifier"])
				id <- paste(id, collapse="#") 
				dtab[i,"FusionID"] <- id
				dtab[i,"FusionSC"] <- sup 
			 }
			 if(sum(grepl("fusion2", colnames(m)))) {
				sup <- unique(m[, "superclass_name"]) 
				sup <- paste(sup, collapse="#") 
				id <- unique(m[, "Identifier"])
				id <- paste(id, collapse="#") 
				dtab[i,"ConsensusID"] <- id
				dtab[i,"ConsensusSC"] <- sup 
			 }
			 mlist[[i]] <- m
		} else {
			 hc <- hclust(dist(mds), method="ward.D")

			 dtree <- cutreeDynamic(dendro = hc, cutHeight = NULL, minClusterSize = 2, method = "hybrid", deepSplit = 2, distM = as.matrix(dist(mds)))
			if(sum(dtree==0)) dtree <- dtree+1

			#cdt <- data.frame(Dim_1=mds[,1], Dim_2=mds[,2], clustid=as.factor(dtree))

			mtmp <- as.matrix(database[match(m[,1], database[,4]), c("superclass_name", "class_name")], nrow=nrow(m))
			m <- cbind(m[,1], matrix('', ncol=1, nrow=nrow(m)),  mtmp,  m[,-1] )
			m <- m[order(dtree), ]
			rownames(m) <- NULL 
			colnames(m)[1:2] <- c("Identifier", "MCSS")

			 m[,2] <- dtree[order(dtree)] 
			 if(sum(grepl("Score", colnames(m)))) {
				psc <- order(as.numeric(m[,"Score"]), decreasing=TRUE)[1]
				sup <- unique(m[which(m[, "MCSS"] == m[psc,"MCSS"]), "superclass_name"]) 
				sup <- sup[sup!=""] 
				sup <- paste(sup, collapse="#") 
				id <- unique(m[which(m[, "MCSS"] == m[psc,"MCSS"]), "Identifier"])
				id <- paste(id, collapse="#") 
				dtab[i,"MetFragID"] <- id
				dtab[i,"MetFragSC"] <- sup 
			 }
			 if(any(colnames(m)=="fusion")) {
				psc <- order(as.numeric(m[,"fusion"]), decreasing=TRUE)[1]
				sup <- unique(m[which(m[, "MCSS"] == m[psc,"MCSS"]), "superclass_name"]) 
				sup <- sup[sup!=""] 
				sup <- paste(sup, collapse="#") 
				id <- unique(m[which(m[, "MCSS"] == m[psc,"MCSS"]), "Identifier"])
				id <- paste(id, collapse="#") 
				dtab[i,"FusionID"] <- id
				dtab[i,"FusionSC"] <- sup 
			 }
			 if(sum(grepl("fusion2", colnames(m)))) {
				psc <- order(as.numeric(m[,"fusion2"]), decreasing=TRUE)[1]
				sup <- unique(m[which(m[, "MCSS"] == m[psc,"MCSS"]), "superclass_name"]) 
				sup <- sup[sup!=""] 
				sup <- paste(sup, collapse="#") 
				id <- unique(m[which(m[, "MCSS"] == m[psc,"MCSS"]), "Identifier"])
				id <- paste(id, collapse="#") 
				dtab[i,"ConsensusID"] <- id
				dtab[i,"ConsensusSC"] <- sup 
			 }
			mlist[[i]] <- m
		}
	  } else {
		 mtmp <- as.matrix(database[match(m[,1], database[,4]), c("superclass_name", "class_name")], nrow=nrow(m))
		 m <- cbind(matrix(m[,1], nrow=nrow(m)), matrix('', ncol=1, nrow=nrow(m)),  mtmp,  matrix(m[,-1], nrow=nrow(m)) )
		 if(sum(grepl("Score", colnames(m)))) {
			sup <- unique(m[, "superclass_name"]) 
			sup <- paste(sup, collapse="#") 
			id <- unique(m[, "Identifier"])
			id <- paste(id, collapse="#") 
			dtab[i,"MetFragID"] <- id
			dtab[i,"MetFragSC"] <- sup 
		 }
		 if(any(colnames(m)=="fusion")) {
			sup <- unique(m[, "superclass_name"]) 
			sup <- paste(sup, collapse="#") 
			id <- unique(m[, "Identifier"])
			id <- paste(id, collapse="#") 
			dtab[i,"FusionID"] <- id
			dtab[i,"FusionSC"] <- sup 
		 }
		 if(sum(grepl("fusion2", colnames(m)))) {
			sup <- unique(m[, "superclass_name"]) 
			sup <- paste(sup, collapse="#") 
			id <- unique(m[, "Identifier"])
			id <- paste(id, collapse="#") 
			dtab[i,"ConsensusID"] <- id
			dtab[i,"ConsensusSC"] <- sup 
		 }
		 mlist[[i]] <- m
	  }
	}
        write.table(dtab, 'metfrag_out/final_summary.txt', sep="\t", row.names=FALSE, quote=FALSE)
        save(mlist, file="metfrag_out/mlist.rda")

