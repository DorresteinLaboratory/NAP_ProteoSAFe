#!/PATHCONDAENV/envs/nap/bin/Rscript

suppressPackageStartupMessages(library(optparse))

 option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  # If a file is given, take id from param file
  make_option(c("-d", "--dir"), type="character", default="",
              help="directory containing fusion scoring info")

 )     
  opt <- parse_args(OptionParser(option_list=option_list))

	load_obj <- function(f)	{
	    env <- new.env()
	    nm <- load(f, env)[1]
	    env[[nm]]
	}

	print(paste("Here is the input:", opt$dir))

	dir <- strsplit(opt$dir, "/")[[1]][1]

	fls <- dir(dir, full.name=TRUE)
	fls <- fls[grep("\\.rda", fls)]

	targL <- list()

	targL[as.numeric(gsub("\\D", "", fls))] <- lapply(fls, load_obj)
	
	if(grepl("consensus", opt$dir)) {
		print(paste("Length of targL:", length(targL)))
                opt$out <- "merge_consensus/consensus.rda"
		targL2 <- targL
		save(targL2, file=opt$out)
	} else {
		print(paste("Length of targL:", length(targL)))
                opt$out <- "merge_fusion/fusion.rda"
		save(targL, file=opt$out)
		
		print(paste("This is the folder after saving:", dir(opt$dir)))
		write.table("All done.", "fusion_res/jobid.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
	}


