#!/PATHCONDAENV/envs/nap/bin/Rscript

suppressPackageStartupMessages(library(optparse))
source("/PATHTOOLSFOLDER/nap_ccms2/Snap/code/get.metfrag.R")

 option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  # If a file is given, take id from param file
  make_option(c("-f", "--file"), type="character", default="",
              help="rda binary file with all variables"),
  make_option(c("-O", "--out"), type="character", default="",
              help="Mock output")

 )     
  opt <- parse_args(OptionParser(option_list=option_list))
  # Load all requirements from previous step
  load(opt$file)
  get.metfrag(query, allspectra, database, type, path, abs.diff = 0.1, ppm, Mdiff, sl_param, IsPositiveIonMode, PrecursorIonMode, idx, sub("\\.psv$", "", opt$out)) 

