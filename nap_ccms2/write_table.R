#!/PATHCONDAENV/envs/nap/bin/Rscript

library(rjson) 

jtab <- fromJSON(file='split_data/config.json')
args <- commandArgs(TRUE)

	  load("split_data/tabgnps.rda")
	  load("lid_res/lid.rda")
	  load("merge_fusion/fusion.rda")
	  load("merge_consensus/consensus.rda")
	  jobid <- readLines("split_data/jobid.txt")

		metfp <- c()
		fusfp <- c()
		fus2fp <- c()
		# This block generates one table for each
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
		  for(k in 1:n) {
				if(sum(colnames(targtemp)=="Score")) {	
					if(as.numeric(targtemp[k, "Score"])==1) {	
						metfp <- c(metfp, paste(targtemp[k,"Identifier"], tabgnps1[i,1], sep="_"))
					}
				}
				if(sum(colnames(targL[[i]])=="fusion")) {	
					if(as.numeric(targL[[i]][k, "fusion"])==1) {	
						fusfp <- c(fusfp, paste(targL[[i]][k,"Identifier"], tabgnps1[i,1], sep="_"))
					}
				}
				if(sum(colnames(targL2[[i]])=="fusion2")) {	
					if(as.numeric(targL2[[i]][k, "fusion2"])==1) {	
						fus2fp <- c(fus2fp, paste(targL2[[i]][k,"Identifier"], tabgnps1[i,1], sep="_"))
					}
				}
			  
			 
		  }
		}

		mdp <- which(duplicated(sub("^.+_", "", metfp)))
		if(length(mdp)) metfp <- metfp[-mdp]
		fdp <- which(duplicated(sub("^.+_", "", fusfp)))
		if(length(fdp)) fusfp <- fusfp[-fdp]
		f2dp <- which(duplicated(sub("^.+_", "", fus2fp)))
		if(length(f2dp)) fus2fp <- fus2fp[-f2dp]

		m <- tabgnps1[,c("cluster.index", "parent.mass", "RTMean", "LibraryID")]

		m <- matrix(apply(m, 2, unlist), ncol=4)
		m <- gsub("\"", "", m)
		m <- cbind(m, matrix("", nrow=nrow(m), ncol=10))
		colnames(m) <- c("cluster.index", "parent.mass", "RTMean", "LibraryID", "MetFragID", "MetFragStruct", "FusionStruct", "ConsensusStruct", "MetFragNode", "FusionNode", "ConsensusNode", "MetFragFrag", "FusionFrag", "ConsensusFrag")
		m[,"MetFragID"] <- "link" 
		m[,9:14] <- m[,1] 

		colnames(m) <- c("cluster.index", "parent.mass", "RTMean", "LibraryID", "MetFragID", "MetFragStruct", "FusionStruct", "ConsensusStruct", "MetFragNode", "FusionNode", "ConsensusNode", "MetFragFrag", "FusionFrag", "ConsensusFrag")

		write.table(m,
			args[7], row.names=FALSE, quote=FALSE, sep="\t") 	
		
		write.table(read.delim('metfrag_out/final_summary.txt'), args[8], row.names=FALSE, quote=FALSE, sep="\t") 
		zip(args[9], c(
				"split_data/tabgnps.rda",
				"split_data/net.rda",
				"split_data/allspectra.rda",
				"merge_fusion/fusion.rda",
				"merge_consensus/consensus.rda"
				),
			flags="-r"
		)

