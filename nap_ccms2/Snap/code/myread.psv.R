myread.psv <- function(fn) {
  mols = try(as.matrix(read.csv(fn, sep="|")), TRUE)
  if(!length(grep("Error", mols))) {
    return(mols)
  } else {
    return(NULL)
  }
}
