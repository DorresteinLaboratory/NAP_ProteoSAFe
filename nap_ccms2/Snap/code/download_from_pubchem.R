

library(rjson) 
#library(httr) 

download_candidates <- function(mass, ppm) {
	mdiff <- (ppm/(10^6))*mass

	mimass <- mass - mdiff
	mxmass <- mass + mdiff
	url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&term=", mimass, "[MIMass]:", mxmass, "[MIMass]&retmode=json&retmax=100000")

	idlist <- fromJSON(file=url)$esearchresult$idlist
	if(!length(idlist)) return(NULL)
	idlist <- split(idlist, ceiling(seq_along(idlist)/50))

	get.properties <- function(idlist) {
		cidlist <- paste(idlist, collapse=",")
		url2 <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cidlist, "/property/MonoisotopicMass,InChI,CanonicalSMILES,InChIKey,MolecularFormula/JSON")
		res <- fromJSON(file=url2)
		do.call(rbind, res$PropertyTable$Properties) 
	}

	fres <- lapply(idlist, get.properties)
	# check exception
	if(is.null(fres)) return(NULL)
	fres <- do.call(rbind, fres) 

	if(nrow(fres)==1) {
		cn <- colnames(fres)	
		fres <- matrix(fres, nrow=1)
		colnames(fres) <- cn
	} else {
		fres <- apply(fres, 2, unlist) 
	}
	
	

	fres <- data.frame(fres, kingdom_name="", superclass_name="", class_name="", subclass_name="") 
	colnames(fres)[1] <- "Identifier" 
	colnames(fres)[3] <- "SMILES" 

	ktab <- sapply(as.character(fres$InChIKey), strsplit, "-") 

	fres$InChIKey2 <- unlist(lapply(ktab, function(x) x[2])) 
	fres$InChIKey1 <- unlist(lapply(ktab, function(x) x[1]))

	fres <- fres[,-5]

	cnames <- c("MonoisotopicMass", "InChI", "SMILES", "Identifier", "InChIKey2", "InChIKey1", "MolecularFormula", "kingdom_name", "superclass_name", "class_name", "subclass_name")

	fres <- fres[, cnames]

#	payload <- list(query_input= fres[1,4]$InChI, query_type='STRUCTURE', label= "Ch1")
#
#	req <- httr::POST("http://classyfire.wishartlab.com/queries.json",
#	  httr::add_headers(
#	     "Content-Type" =  "application/json"
#	  ),
#	  body = payload,
#	  encode = "json"
#	   
#	);
#
#	json <- httr::content(req, as = "text") 
#	qids <- fromJSON(json) 
#
#	sget <- fromJSON(file=paste0("http://classyfire.wishartlab.com/queries/", qids$id, ".json"))
	return(fres)
}
