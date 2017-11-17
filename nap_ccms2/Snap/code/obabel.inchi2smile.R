obabel.inchi2smile <- function(inchi) {
  write.table(inchi, "tmp.inchi", row.names = FALSE, col.names = FALSE, quote = FALSE)
  system("/PATHCONDAENV/envs/nap/bin/obabel tmp.inchi -O tmp.smi")
  res <- sub("\t", "", readLines("tmp.smi"))
  file.remove("tmp.smi")
  file.remove("tmp.inchi")
  return(res)
}
