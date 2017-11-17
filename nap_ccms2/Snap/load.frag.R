#!/PATHCONDAENV/envs/nap/bin/Rscript

suppressPackageStartupMessages(library(optparse))
source("/PATHTOOLSFOLDER/nap_ccms2/Snap/code/myread.psv.R")

	# Ideally should take parameters from command line
        dir <- "fragmenter_res"
        load("split_data/tabgnps.rda")

	fls <- dir(dir, full.name=TRUE)
	fls <- fls[grep("\\.psv", fls)]
	fls <- fls[order(as.numeric(gsub("\\D", "", fls)))]

	lid <- lapply(fls, myread.psv)

	lid2 <- lid

	lid <- as.list(rep(0, nrow(tabgnps1)))
	lid[as.numeric(gsub("\\D", "", fls))] <- lid2
	
	# The workflow manager should create a directory
	save(lid, file="lid_res/lid.rda")


