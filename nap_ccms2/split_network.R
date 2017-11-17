#!/PATHCONDAENV/envs/nap/bin/Rscript


# Set path for java 
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ':/PATHCONDAENV/envs/nap/bin/'))
Sys.setenv(JAVA_HOME = '/PATHCONDAENV/envs/nap/jre')
Sys.setenv(LD_LIBRARY_PATH = '/PATHCONDAENV/envs/nap/jre/lib/amd64:/PATHCONDAENV/envs/nap/jre/lib/amd64/server')

suppressPackageStartupMessages(library(doMC))
library(rcdk) 
library(rjson) 
library(graph)

	# Ideally the parameters are set on command line
	opt <- list()
	opt$file <- "" 
	#opt$charge <- "+" 
	#opt$adduct <- "[M+H]" 
	#opt$database <- "mock" 
	#opt$cnun <- 0
	#opt$ppm <- 15
	#opt$path <- ""
	opt$output <- "split_data"
	
	# Read and parse parameter file
	altargs <- commandArgs(TRUE)
	altargtab <- readLines(altargs[1])
	
	opt$adduct <- sub("<parameter name=.+>(.+)</parameter>", "\\1", altargtab[grep("ADDUCT", altargtab)])
	opt$charge <- sub("<parameter name=.+>(.+)</parameter>", "\\1", altargtab[grep("CHARGE", altargtab)])
	opt$database <- sub("<parameter name=.+>(.+)</parameter>", "\\1", altargtab[grep("DATABASE", altargtab)])
	opt$ppm <- sub("<parameter name=.+>(.+)</parameter>", "\\1", altargtab[grep("PPM", altargtab)])

	opt$job <- sub("<parameter name=.+>(.+)</parameter>", "\\1", altargtab[grep("JOBID", altargtab)])
	opt$cnun <- sub("<parameter name=.+>(.+)</parameter>", "\\1", altargtab[grep("CNUN", altargtab)])
	opt$cosine <- as.numeric(sub("<parameter name=.+>(.+)</parameter>", "\\1", altargtab[grep("COSINE", altargtab)]))
	opt$nfirst <- sub("<parameter name=.+>(.+)</parameter>", "\\1", altargtab[grep("NFIRST", altargtab)])
	opt$prop <- sub("<parameter name=.+>(.+)</parameter>", "\\1", altargtab[grep("PROP", altargtab)])
        opt$class <- sub("<parameter name=.+>(.+)</parameter>", "\\1", altargtab[grep("CCCLASS", altargtab)]) 
        opt$skip <- sub("<parameter name=.+>(.+)</parameter>", "\\1", altargtab[grep("SKIP", altargtab)]) 
        opt$PLOT <- sub("<parameter name=.+>(.+)</parameter>", "\\1", altargtab[grep("PLOT", altargtab)]) 
	opt$nreport <- as.numeric(sub("<parameter name=.+>(.+)</parameter>", "\\1", altargtab[grep("NREPORT", altargtab)]))

        opt$madduct <- sub("<parameter name=.+>(.+)</parameter>", "\\1", altargtab[grep("MULTION", altargtab)]) 

	if(length(opt$prop)) {
		opt$prop <- 1 
	} else {
		opt$prop <- 0 
	}


        if(any(grepl("spec1", altargs))) {
        	opt$udb <- altargs[2]
	} else {
		opt$udb <- "0"
        }

	if(length(opt$skip)) {
		opt$skip <- 1 
	} else {
		opt$skip <- 0 
	}


	PATH <- "/PATHTOOLSFOLDER/nap_ccms2/Snap/"
	if(PATH!=""){
		opt$path <- PATH
	}  

	source(paste0(opt$path, "code/obabel.inchi2smile.R")) 
	source(paste0(opt$path, "code/load.db.R")) 
	source(paste0(opt$path, "code/download_from_pubchem.R")) 

	# Download GNPS network
	fname <- paste(sample(letters, 10), collapse="") 
	if(opt$job!="0") {
		  system(paste("/PATHCONDAENV/envs/nap/bin/python", paste0(opt$path, "code/download_gnps.py -j"), opt$job, "-t Network")) 
		  dir.create(fname)
		  unzip("download.zip", exdir=fname)
		  file.remove("download.zip")
	} else if (opt$file!="" & opt$job=="0") { 
		  dir.create(fname)
		  unzip(dir(opt$file, recursive=TRUE, full.names=TRUE), exdir=fname)
	} else {
		stop("No file or job id")
	}

	# We need a reference table for gnps attribute table, with library ids and structural info
	# Do new version include a header?
	# Network edge list
	net <- read.delim(dir(paste0(fname, "/networkedges/"), full.names=TRUE), header = FALSE, quote="")
	gnps <- read.delim(dir(paste0(fname, "/clusterinfosummarygroup_attributes_withIDs/"), full.names=TRUE), quote="", stringsAsFactors=FALSE)
	
	# Check alteration
	gnps <- gnps[which(gnps[,1] %in% unique(as.vector(as.matrix(net[,1:2])))),]

	if (opt$cnun!=0) { 
		# Retrieve the connected component
		# and disconnect the edges when necessary 
		cat("Searching connected components ...", "\n")
		# If we want to follow the order presented 
		# at web interface, just order by numer of nodes
		conn <- by(net, net[,7], function(x) unique(as.vector(as.matrix(x[,1:2]))))
		cat("Found ", length(conn), "connected components", "\n")
		pclust <- which(unlist(lapply(conn, function(x) sum(x==opt$cnun)))==1)
		comp <- as.numeric(conn[[pclust]]) 
		net <- net[(which(net[,1] %in% comp & net[,2] %in% comp)),] 
	      if(sum(net[,5] > opt$cosine)) {
			net <- net[net[,5]> opt$cosine,]
			gr <- ftM2graphNEL(as.matrix(net[,1:2]))
			conn <- connComp(gr)
			p <- which(unlist(lapply(conn, function(x) sum(x==opt$cnun)))==1)
			comp <- as.numeric(conn[[p]]) 
			net <- net[(which(net[,1] %in% comp & net[,2] %in% comp)),] 
	      }
	      gnps <- gnps[which(gnps[,1] %in% comp),] 
	}
	# The files below have to be extracted from GNPS output  
	# Download the file, extract it and process each component

	# Need GNPS library enriched with struct info  
	# this should help structures on GNPS that do not contain 
	# structure info associated
	tabgnps0 <- read.delim(paste0("/PATHTOOLSFOLDER/nap_ccms2/DB/gnps_network_structural_information_update_only_inchi_smile_casUPDATE08_2016.tsv"), stringsAsFactors=FALSE)

	# Read all spectra
	# Not very efficient way to parse mgf
	mgf.spectra <- dir(fname, full.names=TRUE)
	mgf.spectra <- mgf.spectra[grep("\\.mgf", mgf.spectra)]  

	if(opt$job=="") {
		jfile <- dir(fname, full.names=TRUE)
		jfile <- jfile[grep("\\.xml", jfile)] 
		jobid <- readLines(jfile) 
		jobid <- sub("<parameter name=\"task\">(.+)</parameter>", "\\1", jobid[grep("task", jobid)])
		opt$job <- jobid
	}


        if(any(grepl("spec2", altargs))) {
		sl_param <- readLines(dir("spec2", full.names=TRUE))
	} else {
		sl_param <- readLines(paste0(opt$path, "sl_mock_parameter.txt"))
        }
        

	# Need new file for library ids containing unique library identifier
	specnets <- read.delim(dir(paste0(fname, "/result_specnets_DB/"), full.names=TRUE), quote="", stringsAsFactors=FALSE)
	specnets <- specnets[, c("LibMZ", "Precursor_MZ", "Smiles", "SpectrumID", "MassDiff", "SpecMZ", "MZErrorPPM", "X.Scan.", "INCHI")]

	cids <- unique(as.character(specnets[specnets$SpectrumID %in% as.character(tabgnps0$SpectrumID), "SpectrumID"]))
	# change to match function
	if(length(cids)) {
		cidsinspec <- sapply(cids, function(x) which(specnets$SpectrumID==x)[1])
		cidsintab <- sapply(cids, function(x) which(tabgnps0[,"SpectrumID"]==x)[1])

		tabgnps0$Smiles <- as.character(tabgnps0$Smiles)
		tabgnps0$INCHI <- as.character(tabgnps0$INCHI)
		specnets$INCHI <- as.character(specnets$INCHI)
		specnets$Smiles <- as.character(specnets$Smiles)

		for(i in 1:length(cids)) {
		  if(sum(specnets[cidsinspec[[i]], "Smiles"] == "N/A" & specnets[cidsinspec[[i]], "INCHI"] == "N/A")) {
		    specnets[cidsinspec[[i]], c("INCHI", "Smiles")] <- unique(tabgnps0[cidsintab[[i]],-1])[1, c("INCHI", "Smiles")]
		  }
		}

		smi <- sapply(specnets[which(specnets[,"INCHI"]!="N/A" & specnets[,"Smiles"]=="N/A"), "INCHI"], obabel.inchi2smile)
		smi[unlist(lapply(smi, length))==0] <- ""

		specnets[which(specnets[,"INCHI"]!="N/A" & specnets[,"Smiles"]=="N/A"), "Smiles"] <- unlist(smi)
	}

	tabgnps1tmp <- merge(data.frame(gnps, idx=1:nrow(gnps)), data.frame(specnets[, c("X.Scan.", "SpectrumID", "Smiles", "INCHI")], idx2=1:nrow(specnets)), by.x="cluster.index", by.y="X.Scan.")
	tabgnps1 <- gnps
	tabgnps1$SpectrumID <- ""
	tabgnps1$Smiles <- ""
	tabgnps1$INCHI <- ""

	tabgnps1$SpectrumID[tabgnps1tmp$idx] <- as.character(specnets$SpectrumID[tabgnps1tmp$idx2] )
	tabgnps1$INCHI[tabgnps1tmp$idx] <- as.character(specnets$INCHI[tabgnps1tmp$idx2] )

	# Make sure the smiles will not break the code bellow
	smitmp <- as.character(specnets$Smiles[tabgnps1tmp$idx2] )
	smiidx <- sapply(as.character(specnets$Smiles[tabgnps1tmp$idx2] ), function(x) try(get.mol2formula(parse.smiles(x)[[1]]))) 
	if(sum(grepl("Error", smiidx))) smitmp[grep("Error", smiidx)] <- "" 
	tabgnps1$Smiles[tabgnps1tmp$idx] <- smitmp 

	# MetFrag parameter template
	allspectra <- readLines(mgf.spectra)

	# Translate the adduct annotation from GNPS user 
	# to MetFrag adduct form
	# Change the GNPS forms to adapt to MetFrag form
	# Should import adduct prediction in the future
	adducts <-  c("M+H",        "[M+H]",      "[M-H]",        "[M+NH4]",      "M+2H",       "M+Na",       "[M+K]",        "M",          "M+H-H2O",    "M-H2O+H",    "[M+Na]",     "[M+2H]",     "[M+H]+",     "579.171",   
		    "[M]+*",      "[M]*+",      "[M]",       "[M+Na]+",    "[M-H2O+H]+", "[2M+H]+",    "[M+2H]2+",   "[M+H-H2O]+", "[M+NH4]+", "[M+Cl]", "[M+K]")
	adducts2  <- c(1,1,-1, 18,1,23,39,0,1,1,23,1,1,1,0,0,0,23,1,1,1,1,18, 35 )

	tabgnps <- gnps[, c("cluster.index", "precursor.mass") ]
	
	# Load file containing adduct mass
	# FIX M ADDUCT
	if(opt$charge=="Positive") {
	 if(length(opt$madduct)) {
		  madduct <- strsplit(opt$madduct, ",")[[1]] 
		  adductp <- read.csv(paste0("/PATHTOOLSFOLDER/nap_ccms2/DB/ESI-MS-adducts-calculator-Vsept2012_pos.csv"))
		  adductp <- adductp[adductp[,"Charge"]=="1+",]
		  adductp <- adductp[-which(duplicated(adductp[,"Mass"])),] 
		  adductp[,1] <- gsub("\\s", "", as.character(adductp[,1])) 

		  Mdiff <- c() 

		  for(ad in 1:length(madduct)) {
			  Adduct <- madduct[ad] 
			  Mdiff <- c(Mdiff, adductp[which(adductp[,1]==gsub("\\[|\\]", "", Adduct)), "Mass"])
		  }
		  Mdiff <- paste(Mdiff, collapse=";") 

		  tabgnps <- data.frame(tabgnps, IonMode="Positive", Adduct, Mdiff ) 
	 } else{
	  	  madduct <- opt$madduct
		  Adduct <- opt$adduct
		  adductp <- read.csv(paste0("/PATHTOOLSFOLDER/nap_ccms2/DB/ESI-MS-adducts-calculator-Vsept2012_pos.csv"))
		  adductp <- adductp[adductp[,"Charge"]=="1+",]
		  adductp <- adductp[-which(duplicated(adductp[,"Mass"])),] 
		  adductp[,1] <- gsub("\\s", "", as.character(adductp[,1])) 

		  Mdiff <- adductp[which(adductp[,1]==gsub("\\[|\\]", "", Adduct)), "Mass"] 

		  tabgnps <- data.frame(tabgnps, IonMode="Positive", Adduct, Mdiff ) 
	 }
	} else {
	 if(length(opt$madduct)) {
	  	  madduct <- strsplit(opt$madduct, ",")[[1]] 
		  adductn <- read.csv(paste0("/PATHTOOLSFOLDER/nap_ccms2/DB/ESI-MS-adducts-calculator-Vsept2012_neg.csv"))
		  adductn <- adductn[adductn[,"Charge"]=="1-",]
		  adductn <- adductn[-which(duplicated(adductn[,"Mass"])),] 
		  adductn[,1] <- gsub("\\s", "", as.character(adductn[,1])) 

		  Mdiff <- c() 

		  for(ad in 1:length(madduct)) {
			  Adduct <- madduct[ad] 
			  Mdiff <- c(Mdiff, adductn[which(adductn[,1]==gsub("\\[|\\]", "", Adduct)), "Mass"])
		  }
		  Mdiff <- paste(Mdiff, collapse=";") 

		  tabgnps <- data.frame(tabgnps, IonMode="Negative", Adduct, Mdiff ) 
	
	 } else{
	  	  madduct <- opt$madduct
		  Adduct <- opt$adduct
		  adductn <- read.csv(paste0("/PATHTOOLSFOLDER/nap_ccms2/DB/ESI-MS-adducts-calculator-Vsept2012_neg.csv"))
		  adductn <- adductn[adductn[,"Charge"]=="1-",]
		  adductn <- adductn[-which(duplicated(adductn[,"Mass"])),] 
		  adductn[,1] <- gsub("\\s", "", as.character(adductn[,1])) 

		  Mdiff <- adductn[which(adductn[,1]==gsub("\\[|\\]", "", Adduct)), "Mass"] 

		  tabgnps <- data.frame(tabgnps, IonMode="Negative", Adduct, Mdiff ) 
	 }
	}

	# Load database combination 
#	if(opt$database=="PubChem" | opt$database=="PubChem[Only]") {
#		  TYPE <- "PubChem"
#		  database <- data.frame(MonoisotopicMass=0, InChI="", SMILES="", Identifier="", InChIKey2="", 
#				InChIKey1="", MolecularFormula="", kingdom_name="", superclass_name="", class_name="", subclass_name="")
#	} else {
	 ppm <- as.numeric(opt$ppm)
		  TYPE <- "local"
		  if(opt$database=="none" & opt$udb!="0") {
			database = fread(opt$udb)
			write.table(database, "split_data/udb.txt", sep="\t", row.names=FALSE)

		  } else if(opt$database!="none" & opt$udb!="0") {
			pubchem <- ""
			databasep <- NULL
			dbs <- strsplit(opt$database, ",")[[1]]
			if(any(grepl("PubChem", dbs))) {
				dbs <- dbs[-grep("PubChem", dbs)]
				pubchem <- "PubChem"
			}
			

			if(pubchem=="PubChem") {
		  		databaseL <- list()
					
				for(p in 1:nrow(tabgnps)) {
				  if(tabgnps[p,"IonMode"]=="Positive") {
				        databaseL[[p]] <- download_candidates(as.matrix(tabgnps[p, "precursor.mass"]) - tabgnps[p,"Mdiff"], ppm)
				  } else {
				      if(tabgnps[p,"Adduct"]!="[M+Cl]") {
					databaseL[[p]] <- download_candidates(as.matrix(tabgnps[p, "precursor.mass"]) - tabgnps[p,"Mdiff"], ppm)
				      } else {
					databaseL[[p]] <- download_candidates(as.matrix(tabgnps[p, "precursor.mass"]) - tabgnps[p,"Mdiff"], ppm)
				      }
				  }
				}
				databasep <- do.call(rbind, databaseL)
		  	}

			database2 = fread(opt$udb)
			if(!is.null(databasep)) {
				database2 = rbind(databasep, database2)
			}
			write.table(database2, "split_data/udb.txt", sep="\t", row.names=FALSE)

			if(length(dbs)) {
				database <- lapply(dbs, function(x) load.db("/PATHTOOLSFOLDER/nap_ccms2/", x))  
				database <- do.call(rbind, database) 
				database = rbind(database, database2)
			} else{
				database = database2
			}

		  } else {
			pubchem <- ""
			databasep <- NULL
			dbs <- strsplit(opt$database, ",")[[1]]
			if(any(grepl("PubChem", dbs))) {
				dbs <- dbs[-grep("PubChem", dbs)]
				pubchem <- "PubChem"
			}
			

			if(pubchem=="PubChem") {
		  		databaseL <- list()
					
				for(p in 1:nrow(tabgnps)) {
				  if(tabgnps[p,"IonMode"]=="Positive") {
				        databaseL[[p]] <- download_candidates(as.matrix(tabgnps[p, "precursor.mass"]) - tabgnps[p,"Mdiff"], ppm)
				  } else {
				      if(tabgnps[p,"Adduct"]!="[M+Cl]") {
					databaseL[[p]] <- download_candidates(as.matrix(tabgnps[p, "precursor.mass"]) - tabgnps[p,"Mdiff"], ppm)
				      } else {
					databaseL[[p]] <- download_candidates(as.matrix(tabgnps[p, "precursor.mass"]) - tabgnps[p,"Mdiff"], ppm)
				      }
				  }
				}
				databasep <- do.call(rbind, databaseL)
		  	}
			if(!is.null(databasep)) {
				opt$udb <- 1 
				write.table(databasep, "split_data/udb.txt", sep="\t", row.names=FALSE)
			} else {
		  		databasep <- data.frame(MonoisotopicMass=0, InChI="", SMILES="", Identifier="", InChIKey2="", 
					InChIKey1="", MolecularFormula="", kingdom_name="", superclass_name="", class_name="", subclass_name="")
			}
			if(length(dbs)) {
				database <- lapply(dbs, function(x) load.db("/PATHTOOLSFOLDER/nap_ccms2/", x))  
				database <- do.call(rbind, database) 
				database = rbind(database, databasep)
			} else{
				database <- databasep
			}

		  }
	 	database <- as.data.frame(database) 
	 	database[, 1] <- as.numeric(as.matrix(database[,1])) 
		print(nrow(database))
	 #}

         print(altargs)
         if(any(grepl("spec2", altargs))) {
		mpar <- readLines(altargs[3])
		write.table(mpar, "split_data/sl_mock_parameter.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
	 }

	 system(paste("rm -rf", fname)) 
	
	 # Set variables to be used in the next steps
	 path <- paste0(opt$path, "MetFrag2.3-CL.jar")
	 type <- TYPE
	 output <- strsplit(opt$output, "/")[[1]][1] 
	 database2 <- database
	 allspectra2 <- allspectra
	 rm(database) 
	 rm(allspectra)
	 if(length(opt$class)) {
		class <- opt$class
	 } else {
		class <- NULL
	 }
	 skip <- opt$skip

	dir.create(output)
	setwd(output)
	json <- toJSON( opt )
	write.table(json, "config.json", row.names=FALSE, col.names=FALSE, quote=FALSE)
	print(madduct)
	print(tabgnps[1,])

	# Create Rdata with all required variables to MetFrag step
	do.line <- function(tabgnps, i, database2, ppm, allspectra2, sl_param, path, adducts, adducts2, type, class, skip, madduct) { 
		 if(tabgnps[i,"IonMode"]=="Positive") IsPositiveIonMode <- TRUE else IsPositiveIonMode <- FALSE
		 if(length(which(adducts==tabgnps[i,"Adduct"]))) PrecursorIonMode <- adducts2[which(adducts==tabgnps[i,"Adduct"])] else PrecursorIonMode <- 1
		 query <- tabgnps[i,]
		 Mdiff <- tabgnps[i,"Mdiff"]
		 idx = i

		 if(skip==0) {
			   if(!is.null(class) & !is.null(database2)) {
				cname <- strsplit(class, ":")[[1]][1] 
				cls <- strsplit(class, ":")[[1]][2] 
		      		id <- grepl(cls, database2[,cname])
		      		if(sum(id)) {
					database2 <- database2[id,]
				}
		  	   }	
			   vhmdb <- as.numeric(database2[, 1])
			   if (ppm != 0) {
				 # If query has the exact mass, should be done before
				 if(IsPositiveIonMode) {
	 				if(length(madduct)) {
						Mdiff <- as.numeric(strsplit(as.character(Mdiff), ";")[[1]])
				      		hid <- sapply(Mdiff, function(x) which((10^6 * (abs(vhmdb - (as.matrix(query[2]) - x))/vhmdb)) < ppm))
						hid <- unlist(hid)
					} else{ 
				      		hid <- which((10^6 * (abs(vhmdb - (as.matrix(query[2]) - Mdiff))/vhmdb)) < ppm)
					}
				  } else {
	 				if(length(madduct)) {
						Mdiff <- as.numeric(strsplit(as.character(Mdiff), ";")[[1]])
				      		hid <- sapply(Mdiff, function(x) which((10^6 * (abs(vhmdb - (as.matrix(query[2]) - x))/vhmdb)) < ppm))
						hid <- unlist(hid)
					} else{
					      if(tabgnps[i,"Adduct"]!="[M+Cl]") {
						hid <- which((10^6 * (abs(vhmdb - (as.matrix(query[2]) - Mdiff))/vhmdb)) < ppm)
					      } else {
						hid <- which((10^6 * (abs(vhmdb - (as.matrix(query[2]) - Mdiff))/vhmdb)) < ppm)
					      }
					}
				  }
			    }
			    if(length(hid)) {
				database <- database2[hid, ]
			    } else {
				database <- NULL 
			    }
		  } else {
				database <- database2
		  }
		  if(!is.null(class) & !is.null(database) & skip!=0) {
			cname <- strsplit(class, ":")[[1]][1] 
			cls <- strsplit(class, ":")[[1]][2] 
		        id <- grepl(cls, database[,cname])
		        if(sum(id)) {
				database <- database[id,]
				if(sum(database[,1]>1000)) database <- database[database[,1]<1000,]
			} else {
				database <- NULL 
			}
		  }
		  begin <- which(paste0("SCANS=", as.matrix(query[1])) == allspectra2)
		  end <- grep("END IONS", allspectra2[begin:length(allspectra2)])[1]
		  allspectra <- allspectra2[begin:(begin + end)]   

		 save(list= c("query", "allspectra", "database", "type", "path", "ppm", "Mdiff", "sl_param", "IsPositiveIonMode", "PrecursorIonMode", "idx"), file=paste0("line", i, ".rda") )
	}

#	if(type!="PubChem") {
		#no_cores <- detectCores()
		#registerDoMC(ceiling(no_cores/2))
		#foreach(x=1:nrow(tabgnps)) %dopar% {
		for(x in 1:nrow(tabgnps)) {
			do.line(tabgnps, x, database2, ppm, allspectra2, sl_param, path, adducts, adducts2, type, class, skip, madduct)
		}
#	} else {
#		lerror <- list()
#		for(x in 1:nrow(tabgnps)) {
#			lerror[[x]] <- try(do.line(tabgnps, x, database2, ppm, allspectra2, sl_param, path, adducts, adducts2, type, class, skip, madduct), TRUE)
#		}
#		if(any(grepl("Error", lerror))) {
#			verror <- grep("Error", lerror)	
#			for(x in verror) {
#				try(do.line(tabgnps, x, database2, ppm, allspectra2, sl_param, path, adducts, adducts2, type, class, skip, madduct), TRUE)
#			}
#		}
#	}

	fls <- dir(".")
	fls <- fls[grep("line", fls)]
	#fls <- fls[order(as.numeric(gsub("\\D", "", fls)))]

	chunksize = 15
	if(length(fls)/15 > 50) {
		chunksize =  ceiling(length(fls)/50)
	}

	sz <- system("ls -l .", intern=TRUE)
	sztab <- t(sapply(sz[-1], function(x) strsplit(x, " ")[[1]][strsplit(x, " ")[[1]]!=""]))
	rownames(sztab) <- NULL
	sztab <- sztab[grep("line", sztab[,9]),]
	sztab <- data.frame(sztab, parent.mass=tabgnps1[as.numeric(gsub("\\D", "", sztab[,9])), "parent.mass"]) 
	sztab[,5] <- as.numeric(as.matrix(sztab[,5])) 
	# can also be multiplied by a value proportional to size
	if(sum(sztab[,"parent.mass"]> 1000))  sztab[sztab[,"parent.mass"]>1000, 5] <- sztab[sztab[,"parent.mass"]>1000, 5]*100 

	# https://stackoverflow.com/questions/10467579/solving-task-scheduling-or-bin-packing-optimizations-in-r
	assign.job <- function(machines, job) {
	    which.machines <- which.min(lapply(machines, sum))
	    machines[[which.machines]] <- c(machines[[which.machines]], job)
	    machines
	}

	allocate <- function(num.machines, job.times) {
	    machines <- lapply(1:num.machines, function(...) c())
	    Reduce(assign.job,
		   sort(job.times, decreasing=TRUE),
		   machines)
	}

	if(length(fls)>chunksize){
		nchunks <- ceiling(length(fls)/chunksize)
		chunks <- allocate(nchunks, sztab[,5])
	} else {
		# CHECK
		chunks <- list() 
		chunks[[1]] <- sztab[,5]
	}
	for(ch in 1:length(chunks)) {
		mtch <- match(chunks[[ch]], sztab[,5])
		chunks[[ch]] <- as.character(sztab[mtch,9])
		sztab <- sztab[-mtch,]
	}	


	# Save the processed edge list, attribute table and fragmentation spectra for next steps
	save(file=paste0("tabgnps.rda"), list=c("tabgnps1"))
	save(file=paste0("net.rda"), list=c("net"))
	save(file=paste0("allspectra.rda"), list=c("allspectra2"))
	write.table(opt$job, "jobid.txt", row.names=FALSE, col.names=FALSE, quote=FALSE) 
	setwd("..")

	#idxvec <- c(seq(1, length(fls), chunksize), length(fls)+1)
        #for(index in 1:(length(idxvec)-1)) {
        for(index in 1:length(chunks)) {
		#to_save <- fls[idxvec[index]:(idxvec[index+1]-1)]
		write.table(chunks[[index]], paste0("chunk_data/chunk", index, ".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
	}
	t.stop <- as.character(Sys.time())


	cat("Finishing Split at", t.stop, "\n")


