require("data.table")
load.db <- function(path, type) {
       switch(type,
              HMDB = fread(paste0(path, "DB/HMDBCLASS.txt")),
              GNPS = fread(paste0(path, "DB/gnpsCLASS.txt")),
              DNP = fread(paste0(path, "DB/dnpCLASS.txt")),
              CHEBI = fread(paste0(path, "DB/chebiCLASS.txt")),
              SUPNAT = fread(paste0(path, "DB/supernaturalCLASS.txt")),
              mock = fread(paste0(path, "DB/mockCLASS.txt"))
	)
     }

