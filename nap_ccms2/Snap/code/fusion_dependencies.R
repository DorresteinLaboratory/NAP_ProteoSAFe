get.similarity <- function(ind, v1, vidx, tabgnps1) {
  # Check if the associated node has structural information (ID) 
  if(tabgnps1[tabgnps1[,1]==vidx,"Smiles"]!="N/A" & tabgnps1[tabgnps1[,1]==vidx,"Smiles"]!=" " & tabgnps1[tabgnps1[,1]==vidx,"Smiles"]!="") {
    
    # make targ a parameter 
    mols <- c(as.character(tabgnps1[which(tabgnps1[,1]==vidx), "Smiles"]), targ[,"SMILES"])
    # Calculate structural similarity
    tn <- get.tanimoto(mols)
    # attach structural similarity to spectral similarity
    # on column 5 of network matrix
    targ <- cbind(unlist(tn), 
                  v1[which((v1[,1]==tabgnps1[ind,1] & v1[,2]==vidx) | (v1[,2]== tabgnps1[ind,1] & v1[,1]==vidx)),5])
    colnames(targ) <- c("tn", "cosine")
    gc()
    return(targ)
  } else {
    return(NULL)
  }

}

require(rcdk)
get.tanimoto <- function(mols) {
  mols <- parse.smiles(mols)
  invisible(sapply(mols, do.typing))
  invisible(sapply(mols, do.aromaticity))
  fps <- lapply(mols, get.fingerprint, type = "extended")

  nu <- which(unlist(lapply(fps, is.null)))
  if(length(nu)) if(sum(nu==1)) return(as.list(rep(1, (length(fps)-1))) ) 
  if(length(nu)) {
    tn <- as.list(rep(0, length(fps)))
    fps <- fps[-nu]
    tn[-c(1, nu)] <- lapply(fps[-1], function(x) fingerprint::distance(fps[[1]], x, method = "tanimoto"))
    tn <- tn[-1]
  } else {
    tn <- lapply(fps[-1], function(x) fingerprint::distance(fps[[1]], x, method = "tanimoto"))
  }
  tn
}

sc <- function(metfrag, alpha=0.3,  cosine, tn) { 
  alpha*metfrag + (1-alpha)*sum(sigmoid(cosine*tn)) 
}

sigmoid <- function(x, a = 9, b = 0.6) 1/(1+exp(-a*(x-b)))
