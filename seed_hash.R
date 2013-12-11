# Deterministic way to set the seed with arbitrary args

if (!exists('seed.hash')) {

library(digest)

seed.hash <- function (...) {
  dd <- digest(paste(SEED.BASE, ..., sep=":"))
  hh <- paste("0x", substr(dd, 10, 15), sep="")
  as.integer(hh)
}

}
