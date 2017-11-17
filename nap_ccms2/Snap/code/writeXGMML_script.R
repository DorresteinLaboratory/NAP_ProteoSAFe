writeXGMML <- function(stabgnps, net, out) {
	g.new <- '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
	<graph label="NAP" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cy="http://www.cytoscape.org" xmlns="http://www.cs.rpi.edu/XGMML"  directed="1">
	  <att name="documentVersion" value="1.1"/>
	  <att name="networkMetadata">
	    <rdf:RDF>
	      <rdf:Description rdf:about="http://www.cytoscape.org/">
		<dc:type>Protein-Protein Interaction</dc:type>
		<dc:description>N/A</dc:description>
		<dc:identifier>N/A</dc:identifier>
		<dc:date>2016-02-08 11:57:06</dc:date>
		<dc:title>structWindow</dc:title>
		<dc:source>http://www.cytoscape.org/</dc:source>
		<dc:format>Cytoscape-XGMML</dc:format>
	      </rdf:Description>
	    </rdf:RDF>
	  </att>
	  <att type="string" name="backgroundColor" value="#ccccff"/>
	  <att type="real" name="GRAPH_VIEW_ZOOM" value="1.0"/>
	  <att type="real" name="GRAPH_VIEW_CENTER_X" value="0.0"/>
	  <att type="real" name="GRAPH_VIEW_CENTER_Y" value="0.0"/>
	  <att type="boolean" name="NODE_SIZE_LOCKED" value="true"/>
	</graph>'

	# other type: "ELLIPSE"
	s.new <- '<node label="NAME" id="ID">
	    <att type="string" name="canonicalName" value="NAME"/>
	    <graphics type="ROUND_RECTANGLE" outline="#FF0033" fill="#89D0F5" width="3.0" h="50.0" w="50.0">
	    <att name="NODE_TRANSPARENCY" value="0" type="string" cy:type="String"/>
	    <att name="NODE_LABEL" value="" type="string" cy:type="String"/>
	    <att name="NODE_LABEL_TRANSPARENCY" value="0" type="string" cy:type="String"/>
	    <att name="NODE_LABEL_FONT_SIZE" value="1" type="string" cy:type="String"/>
	    <att name="NODE_BORDER_STROKE" value="SOLID" type="string" cy:type="String"/>
	    <att name="NODE_BORDER_TRANSPARENCY" value="255" type="string" cy:type="String"/>
	    <att name="COMPOUND_NODE_SHAPE" value="ROUND_RECTANGLE" type="string" cy:type="String"/>
	    </graphics>
	  </node>'

	e.new <- '  <edge label="NAME" source="IDX1" target="IDX2">
	    <att type="string" name="canonicalName" value="NAME"/>
	    <graphics fill="#848484" width="2.0">
	      <att name="EDGE_LABEL" value="" type="string" cy:type="String"/>
	      <att name="EDGE_LABEL_COLOR" value="#000000" type="string" cy:type="String"/>
	      <att name="EDGE_LABEL_FONT_SIZE" value="10" type="string" cy:type="String"/>
	    </graphics>
	  </edge>'

	n.newback <- '<node label="NAME" id="ID">
	    <att type="string" name="canonicalName" value="NAME"/>
	</node>'

	n.new <- '<node label="NAME" id="ID">
	    <att type="string" name="canonicalName" value="NAME"/>
	  </node>'

	stabgnps$LibraryID <- as.character(stabgnps$LibraryID) 
	stabgnps$LibraryID <- gsub("\"", "", as.character(stabgnps$LibraryID))
	stabgnps$LibraryID <- gsub("<|>", "", as.character(stabgnps$LibraryID))

	attr.mt.final <- as.matrix(stabgnps)
	attr.mt.final <- gsub("\\s", "", attr.mt.final)
	attr.mt.final <- gsub("\"", "", attr.mt.final) 
	f.test.results5 <- net
	  
	g.new <- strsplit(g.new, "\n")[[1]]
	s.new <- strsplit(s.new, "\n")[[1]]
	e.new <- strsplit(e.new, "\n")[[1]]
	n.new <- strsplit(n.new, "\n")[[1]]
	  
	attr.mt.final[,"ProteoSAFeClusterLink"] <- gsub("&", "&amp;", attr.mt.final[,"ProteoSAFeClusterLink"])
	 
	# add color to metfrag/fusion here 
	if(sum(grepl("LibraryID", colnames(attr.mt.final)))) {
		ids <-which(attr.mt.final[,"LibraryID"] !="")
		attr.mt.final <- cbind(attr.mt.final, "#D3D3D3")
		colnames(attr.mt.final)[ncol(attr.mt.final)] <- "fill"
		attr.mt.final[ids,"fill"] <- "#228B22"
	} else {
	    attr.mt.final <- cbind(attr.mt.final, "#D3D3D3")
	    colnames(attr.mt.final)[ncol(attr.mt.final)] <- "fill"
	}
	  
	if(sum(grepl("Smiles", colnames(attr.mt.final), ignore.case = TRUE))) {
		attr.n <- colnames(attr.mt.final)[-c(1, 2, grep("INCHI", colnames(attr.mt.final)), tail(1:ncol(attr.mt.final), 1))]
		attr.n <- colnames(attr.mt.final)[-c(1, 2,  tail(1:ncol(attr.mt.final), 1))]
	} else{
		attr.n <- colnames(attr.mt.final)[-c(1, 2, tail(1:ncol(attr.mt.final), 1))]
	}
	  
        # Need three node types
	# 1 - unknown - circle, solid, parent mass name
	# 2 - library match - square, transparent, green border
	# 3 - in silico only - square, transparent, blue border
	# can use the fact that structures are attributes to add as many as want and print sequentially
	get.att <- function(att.name, v1) {
		tmp <- sub("canonicalName",  att.name, n.new[2])
		sub("NAME", v1[att.name], tmp)
	}
	
	get.node <- function(v1,  attr.n ) {
		if(v1["LibraryID"]!="N/A") {
			n.tmp <- s.new
			n.tmp[1] <- sub("NAME", v1["cluster.index"], n.tmp[1])
			n.tmp[1] <- sub("ID", paste0("-", v1["cluster.index"]), n.tmp[1])
			n.tmp[2] <- sub("NAME", v1["parent.mass"], n.tmp[2])
			n.tmp[3] <- sub("#FF0033", "#008000", n.tmp[3])
		
			vatt <- sapply(attr.n, get.att, v1)
			names(vatt) <- NULL
			
			nnode <- c(n.tmp[1:2],  vatt, n.tmp[3:length(n.tmp)])
			rownames(nnode) <- NULL
		} else if(v1["LibraryID"]=="N/A" & v1["MetFragID"]=="") {
			# 1 - unknown - circle, solid, parent mass name
			n.tmp <- s.new
			n.tmp[1] <- sub("NAME", v1["cluster.index"], n.tmp[1])
			n.tmp[1] <- sub("ID", paste0("-", v1["cluster.index"]), n.tmp[1])
			n.tmp[2] <- sub("NAME", v1["parent.mass"], n.tmp[2])
			n.tmp[3] <- sub("ROUND_RECTANGLE", "ELLIPSE", n.tmp[3])
			n.tmp[3] <- gsub("50", "35", n.tmp[3])
			
			n.tmp[4] <- sub("0", "255", n.tmp[4])
			n.tmp[5] <- sub("value=\"\"", paste0("value=\"", v1["parent.mass"], "\""), n.tmp[5])
			n.tmp[6] <- sub("0", "255", n.tmp[6])
			n.tmp[7] <- sub("1", "10", n.tmp[7])
			n.tmp[9] <- sub("255", "0", n.tmp[9])
			
			vatt <- sapply(attr.n, get.att, v1)
			names(vatt) <- NULL
			
			nnode <- c(n.tmp[1:2],  vatt, n.tmp[3:length(n.tmp)])
			rownames(nnode) <- NULL
	     	} else if(v1["LibraryID"]=="N/A" & v1["MetFragID"]!="") {
			# 3 - in silico only - square, transparent, blue border
			n.tmp <- s.new
			n.tmp[1] <- sub("NAME", v1["cluster.index"], n.tmp[1])
			n.tmp[1] <- sub("ID", paste0("-", v1["cluster.index"]), n.tmp[1])
			n.tmp[2] <- sub("NAME", v1["parent.mass"], n.tmp[2])
			n.tmp[3] <- sub("#FF0033", "#0000ff", n.tmp[3])
			
			vatt <- sapply(attr.n, get.att, v1)
			names(vatt) <- NULL
			
			nnode <- c(n.tmp[1:2],  vatt, n.tmp[3:length(n.tmp)])
			rownames(nnode) <- NULL
	     	}
	    	nnode
	  }
	  
	  all.nodes <- lapply(1:nrow(attr.mt.final), function(x) get.node(attr.mt.final[x,],  attr.n))
	  all.nodes <- do.call(c, all.nodes)
	  
	  get.edge.fomatting <- function(v1) {
		  e.tmp <- e.new
		  tmp <- sub("NAME", paste(v1[1], "(unspecified)", v1[2]), e.tmp[1])
		  tmp <- sub("IDX1", -v1[1], tmp)
		  tmp <- sub("IDX2", -v1[2], tmp)
		    
		  e.tmp[1] <- tmp
		  e.tmp[2] <- sub("NAME", paste(v1[1], "(unspecified)", v1[2]), e.tmp[2])
		  e.tmp[4] <- sub("value=\"\"", paste0("value=\"", v1[3], "\""), e.tmp[4])
		  e.tmp
	  }
	  
	  all.edge <- lapply(1:nrow(f.test.results5), function(x) get.edge.fomatting(f.test.results5[x,]))
	  all.edge <- do.call(c, all.edge)
	  write.table(c(g.new[1:21], all.nodes, all.edge, g.new[22]), out, quote=FALSE, col.names = FALSE, row.names = FALSE)
} 
