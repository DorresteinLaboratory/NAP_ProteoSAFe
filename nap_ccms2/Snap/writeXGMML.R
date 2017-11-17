#!/PATHCONDAENV/envs/nap/bin/Rscript

suppressPackageStartupMessages(library(optparse))
library(rjson)

  source("/PATHTOOLSFOLDER/nap_ccms2/Snap/code/writeXGMML_script.R")

  load("merge_consensus/consensus.rda")
  load("split_data/tabgnps.rda")
  load("split_data/net.rda") 
  load("lid_res/lid.rda")
  load("merge_fusion/fusion.rda")

  jtab <- fromJSON(file='split_data/config.json')

  # Number of candidates structures to report
  nreport <- jtab$nreport 
  stabgnps <- tabgnps1[,c("cluster.index", "parent.mass", "number.of.spectra", "RTMean", "sum.precursor.intensity.", "LibraryID", "SpectrumID", "Smiles", "INCHI", "ProteoSAFeClusterLink")]

stabgnps$MetFragScore <- ""
stabgnps$MetFragSMILES <- ""
stabgnps$MetFragID <- ""

stabgnps$FusionScore <- ""
stabgnps$FusionSMILES <- ""
stabgnps$FusionID <- ""

stabgnps$ConsensusScore <- ""
stabgnps$ConsensusSMILES <- ""
stabgnps$ConsensusID <- ""

for(i in 1:nrow(stabgnps)) {
  if(length(lid[[i]])==1 | is.null(lid[[i]])) next
  if(suppressWarnings(sum(lid[[i]]==0, na.rm = TRUE))!=length(lid[[i]])){
    if(!is.null(nrow(lid[[i]])) & nrow(lid[[i]])> nreport-1 & nrow(lid[[i]])> 1 ) {
      stabgnps[i, "MetFragScore"] <- paste(lid[[i]][1:nreport,"Score"], collapse = ",")  
      stabgnps[i, "MetFragSMILES"] <- paste(lid[[i]][1:nreport,"SMILES"], collapse = ",")  
      stabgnps[i, "MetFragID"] <- paste(lid[[i]][1:nreport,"Identifier"], collapse = ",")  
    } else {
      stabgnps[i, "MetFragScore"] <- paste(lid[[i]][,"Score"], collapse = ",")  
      stabgnps[i, "MetFragSMILES"] <- paste(lid[[i]][,"SMILES"], collapse = ",")  
      stabgnps[i, "MetFragID"] <- paste(lid[[i]][,"Identifier"], collapse = ",")  
    }
  }
  if(i>length(targL2)) next
  if(!is.null(targL2[[i]])){
    if(!is.null(nrow(targL2[[i]])) & nrow(targL2[[i]])> nreport-1 & nrow(targL2[[i]])> 1 ) {
      if(sum(colnames(targL2[[i]])=="fusion")) {
        stabgnps[i, "FusionScore"] <- paste(targL2[[i]][order(as.numeric(targL2[[i]][, "fusion"]), decreasing = TRUE), ][1:nreport,"fusion"], collapse = ",")  
        stabgnps[i, "FusionSMILES"] <- paste(targL2[[i]][order(as.numeric(targL2[[i]][, "fusion"]), decreasing = TRUE), ][1:nreport,"SMILES"], collapse = ",")  
        stabgnps[i, "FusionID"] <- paste(targL2[[i]][order(as.numeric(targL2[[i]][, "fusion"]), decreasing = TRUE), ][1:nreport,"Identifier"], collapse = ",")  
      }
      if(sum(colnames(targL2[[i]])=="fusion2")) {
        stabgnps[i, "ConsensusScore"] <- paste(targL2[[i]][order(as.numeric(targL2[[i]][, "fusion2"]), decreasing = TRUE), ][1:nreport,"fusion2"], collapse = ",")  
        stabgnps[i, "ConsensusSMILES"] <- paste(targL2[[i]][order(as.numeric(targL2[[i]][, "fusion2"]), decreasing = TRUE), ][1:nreport,"SMILES"], collapse = ",")  
        stabgnps[i, "ConsensusID"] <- paste(targL2[[i]][order(as.numeric(targL2[[i]][, "fusion2"]), decreasing = TRUE), ][1:nreport,"Identifier"], collapse = ",")  
      }
    } else {
      if(sum(colnames(targL2[[i]])=="fusion")) {
        if(nrow(targL2[[i]])==1) {
          stabgnps[i, "FusionScore"] <- paste(targL2[[i]][,"fusion"], collapse = ",")  
          stabgnps[i, "FusionSMILES"] <- paste(targL2[[i]][,"SMILES"], collapse = ",")  
          stabgnps[i, "FusionID"] <- paste(targL2[[i]][,"Identifier"], collapse = ",")  
        } else{
          stabgnps[i, "FusionScore"] <- paste(targL2[[i]][order(as.numeric(targL2[[i]][, "fusion"]), decreasing = TRUE), ][,"fusion"], collapse = ",")  
          stabgnps[i, "FusionSMILES"] <- paste(targL2[[i]][order(as.numeric(targL2[[i]][, "fusion"]), decreasing = TRUE), ][,"SMILES"], collapse = ",")  
          stabgnps[i, "FusionID"] <- paste(targL2[[i]][order(as.numeric(targL2[[i]][, "fusion"]), decreasing = TRUE), ][,"Identifier"], collapse = ",")  
        }
      }
      if(sum(colnames(targL2[[i]])=="fusion2")) {
        if(nrow(targL2[[i]])==1) {
          stabgnps[i, "ConsensusScore"] <- paste(targL2[[i]][,"fusion2"], collapse = ",")  
          stabgnps[i, "ConsensusSMILES"] <- paste(targL2[[i]][,"SMILES"], collapse = ",")  
          stabgnps[i, "ConsensusID"] <- paste(targL2[[i]][,"Identifier"], collapse = ",")  
        } else{
          stabgnps[i, "ConsensusScore"] <- paste(targL2[[i]][order(as.numeric(targL2[[i]][, "fusion2"]), decreasing = TRUE), ][,"fusion2"], collapse = ",")  
          stabgnps[i, "ConsensusSMILES"] <- paste(targL2[[i]][order(as.numeric(targL2[[i]][, "fusion2"]), decreasing = TRUE), ][,"SMILES"], collapse = ",")  
          stabgnps[i, "ConsensusID"] <- paste(targL2[[i]][order(as.numeric(targL2[[i]][, "fusion2"]), decreasing = TRUE), ][,"Identifier"], collapse = ",")  
        }
      }
    }
  }
}
  dout  <- "final_out"
  writeXGMML(stabgnps, net, "final_out/structure_graph_alt.xgmml") 
  stabgnps <- gsub("\\s+", "", as.matrix(stabgnps))
  write.table(stabgnps, paste0(dout, "/node_attributes_table.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
