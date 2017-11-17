#!/PATHCONDAENV/envs/nap/bin/Rscript

suppressPackageStartupMessages(library(optparse))

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
  
  # Load functions
  source("/PATHTOOLSFOLDER/nap_ccms2/Snap/code/fusion_dependencies.R")

  load("split_data/tabgnps.rda")
  load("split_data/net.rda")
  load("lid_res/lid.rda")
  opt$out <- sub("fragmenter", "fusion", opt$file) 
  opt$out <- sub("psv$", "rda", opt$out) 

  i <- as.numeric(gsub("\\D", "", opt$file)) 

  if(length(lid[[i]])!=1 & !is.null(lid[[i]])){
    targtemp <- lid[[i]]
  } else {
    tmp <- lid[[i]] 
    save(tmp, file=opt$out)
    q()
  } 
  
  if(is.null(ncol(targtemp))) {
    n <- names(targtemp)
    targtemp <- matrix(targtemp, nrow=1)
    colnames(targtemp) <- n
  }
  
  if(!sum(grepl("Identifier", colnames(targtemp)))){
    	tmp <- lid[[i]] 
    	save(tmp, file=opt$out)
        q()
  }
  targ <- matrix(targtemp[,c("Score", "InChI", "SMILES", "Identifier", "InChIKey1")], ncol=5)
  colnames(targ) <- c("Score", "InChI", "SMILES", "Identifier", "InChIKey1")
  if (targ[1,1]=="0.0" | targ[1,1]=="NA") {
    tmp <- lid[[i]] 
    save(tmp, file=opt$out)
    q()
  } 
  
  # Find compounds connected to unidentified
  v1 <- net[which(net[,1]==tabgnps1[i, 1] | net[,2]==tabgnps1[i, 1]),]
  v2 <- setdiff(unique(as.vector(as.matrix(v1[,1:2]))), tabgnps1[i, 1])
  
  if(length(v2)) {
    # Compute similarity between candidates and structures of all neighbor nodes
    # with spectral library match
    t1 <- lapply(v2, function(x) get.similarity(i, v1, vidx = x, tabgnps1))
    if(sum(unlist(lapply(t1, is.null)))==length(t1)) {
    	save(targ, file=opt$out)
        q()
    }
	
    t1 <- do.call(cbind, t1) 
    colnames(t1) <- paste0(colnames(t1), rep(1:(ncol(t1)/2), each=2))
    targ <- cbind(targ, t1)
   
    if(sum(grepl("cosine", colnames(targ)))) {
      # Calculate fusion score for each candidate
      fusion <- apply(targ, 1, function(x) sc(as.numeric(x[1]), cosine = as.numeric(x[grep("cosine", colnames(targ))]), tn = as.numeric(x[grep("tn", colnames(targ))])))
      # re-scale the score from 0 to 1 
      fusion <- round(fusion/max(fusion), 3)
      targ <- cbind(targ, fusion)
    }
    save(targ, file=opt$out)
    q()
  } else {
    save(targ, file=opt$out)
    q()
  }

