#!/PATHCONDAENV/envs/nap/bin/Rscript

suppressPackageStartupMessages(library(optparse))
library(rjson)

# Each file in fragmenter_res will be used
# as index container

 option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  # If a file is given, take id from param file
  make_option(c("-f", "--file"), type="character", default="",
              help="psv files containing index")
 )     

  opt <- parse_args(OptionParser(option_list=option_list))

  source("/PATHTOOLSFOLDER/nap_ccms2/Snap/code/fusion_dependencies.R")

  load("split_data/tabgnps.rda")
  load("split_data/net.rda")
  load("lid_res/lid.rda")
  load("merge_fusion/fusion.rda")
  plist <- fromJSON(file="split_data/config.json")

  opt$out <- sub("fusion", "consensus", opt$file) 

  i <- as.numeric(gsub("\\D", "", opt$file)) 

  # Use n-first for propagation
  nfirst <- as.numeric(plist$nfirst)
  # Wheter to use the propagation from previous step
  prop <- as.numeric(plist$prop)

  if(i>length(targL)) {
    tmp <- NULL 
    save(tmp, file=opt$out)
    q()
  }
  if(!is.null(targL[[i]])){
    targtemp <- targL[[i]]
  } else {
    tmp <- NULL 
    save(tmp, file=opt$out)
    q()
  } 
  
  if(is.null(ncol(targtemp))) {
    n <- names(targtemp)
    targtemp <- matrix(targtemp, nrow=1)
    colnames(targtemp) <- n
  }
 
  if(!sum(grepl("Identifier", colnames(targtemp)))) {
    tmp <- NULL 
    save(tmp, file=opt$out)
    q()
  }  
  targ <- targtemp
  if (targ[1,1]=="0.0" | targ[1,1]=="NA") {
    tmp <- targ 
    save(tmp, file=opt$out)
    q()
  } 
  
  # Find compounds connected to unidentified
  v1 <- net[which(net[,1]==tabgnps1[i, 1] | net[,2]==tabgnps1[i, 1]),]
  v2 <- setdiff(unique(as.vector(as.matrix(v1[,1:2]))), tabgnps1[i, 1])
  
  if(length(v2)) {
    simL <- list()
    for(j in 1:nrow(v1)) {
      if(which(tabgnps1[, 1] == v1[j, 1]) > length(targL)){
	    tmp <- targ 
	    save(tmp, file=opt$out)
            next
      }
      if(which(tabgnps1[, 1] == v1[j, 2]) > length(targL)){
	    tmp <- targ 
	    save(tmp, file=opt$out)
            next
      }
      if(which(tabgnps1[, 1]==v1[j,1])>length(targL) | is.null(targL[[which(tabgnps1[, 1]==v1[j,1])]]) ) {
	    tmp <- targ 
	    save(tmp, file=opt$out)
            next
      }

      if(which(tabgnps1[, 1]==v1[j,2])>length(targL) | is.null(targL[[which(tabgnps1[, 1]==v1[j,2])]]) ) {
	    tmp <- targ 
	    save(tmp, file=opt$out)
            next
      }

      # Calculate local similarity    
      if(!prop) {
	      m1 <- targL[[which(tabgnps1[, 1]==v1[j,1])]][,"SMILES"]
	      od1 <- 1:length(m1)
	      m2 <- targL[[which(tabgnps1[, 1]==v1[j,2])]][,"SMILES"]
	      od2 <- 1:length(m2)
      } else {
	      # If fusion is used, the candidate list is re-ordered
	      tgtmp1 <- targL[[which(tabgnps1[, 1]==v1[j,1])]]
	      tgtmp2 <- targL[[which(tabgnps1[, 1]==v1[j,2])]]
	      if(is.null(tgtmp1) | length(tgtmp1)==1 | is.null(tgtmp2) | length(tgtmp2)==1){
		    tmp <- targ 
		    save(tmp, file=opt$out)
            	    next
      	       }

	      if(sum(grepl("fusion", colnames(tgtmp1)))){
		od1 <- order(as.numeric(tgtmp1[,"fusion"]), decreasing=TRUE) 
	      	m1 <- tgtmp1[od1,"SMILES"]
	      } else {
		od1 <- 1:nrow(tgtmp1)
	      	m1 <- tgtmp1[,"SMILES"]
	      }
              if(sum(grepl("fusion", colnames(tgtmp2)))){
		od2 <- order(as.numeric(tgtmp2[,"fusion"]), decreasing=TRUE) 
	      	m2 <- tgtmp2[od2,"SMILES"]
	      } else {
	      	m2 <- tgtmp2[,"SMILES"]
		od2 <- 1:nrow(tgtmp2)
	      }
      }

      
      if(v1[j,1]==tabgnps1[i, 1]) {

	      #if(length(m1) > nfirst) m1 <- m1[1:nfirst]
	      if(length(m2) > nfirst) m2 <- m2[1:nfirst]
	      
	      #if(sum(unlist(lapply(m1, length))==0)) m1[unlist(lapply(m1, length))==0] <- ""
	      #if(sum(unlist(lapply(m2, length))==0)) m2[unlist(lapply(m2, length))==0] <- ""
	      
	      if(!length(m1) | !length(m2) | is.null(m1) | is.null(m2)) {
		    tmp <- targ 
		    save(tmp, file=opt$out)
		    next
	      }
	      
	      mols <- lapply(c(m1, m2), function(x) parse.smiles(x)[[1]])
	      fps <- lapply(mols, function(x) try(get.fingerprint(x, type = "extended"), TRUE))
	      if(sum(grepl("Error", fps))){
		    tmp <- targ 
		    save(tmp, file=opt$out)
		    next
	      }
	      fp.sim <- fp.sim.matrix(fps, method = "tanimoto")
	      
	      if(sum(is.nan(fp.sim))) fp.sim[is.nan(fp.sim)] <- 0

        # Take the maximum possible similarity for each edge, hopefully will the better context
	highs <- apply(matrix(fp.sim[1:length(m1), (length(m1)+1):ncol(fp.sim) ], nrow=length(m1)), 1, max)
	if(length(m1) < nfirst) {
		fixpos <- rep(0, length(od1))
		fixpos[od1] <- highs
	} else {
		fixpos <- rep(0, length(od1))
		#fixpos[od1[1:nfirst]] <- highs
		fixpos[od1] <- highs
	}
        simL[[j]] <- cbind(fixpos, v1[j,5])
      } else {

	      if(length(m1) > nfirst) m1 <- m1[1:nfirst]
	      #if(length(m2) > nfirst) m2 <- m2[1:nfirst]
	      
	      #if(sum(unlist(lapply(m1, length))==0)) m1[unlist(lapply(m1, length))==0] <- ""
	      #if(sum(unlist(lapply(m2, length))==0)) m2[unlist(lapply(m2, length))==0] <- ""
	      
	      if(!length(m1) | !length(m2) | is.null(m1) | is.null(m2)) {
		    tmp <- targ 
		    save(tmp, file=opt$out)
		    next
	      }
	      
	      mols <- lapply(c(m1, m2), function(x) parse.smiles(x)[[1]])
	      fps <- lapply(mols, function(x) try(get.fingerprint(x, type = "extended"), TRUE))
	      if(sum(grepl("Error", fps))){
		    tmp <- targ 
		    save(tmp, file=opt$out)
		    next
	      }
	      fp.sim <- fp.sim.matrix(fps, method = "tanimoto")
	      
	      if(sum(is.nan(fp.sim))) fp.sim[is.nan(fp.sim)] <- 0

        

	highs <- apply(matrix(fp.sim[(length(m1)+1):ncol(fp.sim), 1:length(m1)], nrow=length((length(m1)+1):ncol(fp.sim))), 1, max)
	if(length(m2) < nfirst) {
		fixpos <- rep(0, length(od2))
		fixpos[od2] <- highs
	} else {
		fixpos <- rep(0, length(od2))
		#fixpos[od2[1:nfirst]] <- highs
		fixpos[od2] <- highs
	}
        simL[[j]] <- cbind(fixpos, v1[j,5])
      }
    }
    if(!length(simL)) {
	    tmp <- targ 
	    save(tmp, file=opt$out)
	    q() 
    }
    t1 <-  do.call(cbind, simL)
    colnames(t1) <- paste0(rep(c("Sim", "Cos"), ncol(t1)/2), rep(1:(ncol(t1)/2), each=2))
    targ <- cbind(targ, t1)
    
    if(sum(grepl("Cos", colnames(targ)))) {
         # Calculate fusion score
         fusion2 <- apply(targ, 1, function(x) sc(as.numeric(x[1]), cosine = as.numeric(x[grep("Cos", colnames(targ))]), tn = as.numeric(x[grep("Sim", colnames(targ))])))
         # re-scale the score from 0 to 1 
         fusion2 <- round(fusion2/max(fusion2), 3)
         targ <- cbind(targ, fusion2)
    }
    save(targ, file=opt$out)
    q()
  } else {
    save(targ, file=opt$out)
    q()
  }


