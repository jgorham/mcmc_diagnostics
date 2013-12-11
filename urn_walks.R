source('seed_hash.R')
source('mcmc_urn.R')

library(plyr)
library(parallel)

NUM.WORKERS <- 8
#DIM.VALUES <- c(100, 250, 1000, 4000)
DIM.VALUES <- c(4000)
CUT.OFF.MULT <- c(0.5, 1, 2, 4)
SEED.BASE <- '45wrba'

params <- expand.grid(
  d=DIM.VALUES,
  cut.off.mult=CUT.OFF.MULT
)

urn.walks <- adply(params, 1, function (param) {
  print(param)

  d <- param$d
  cut.off.mult = param$cut.off.mult
  theta <- 1
  num.iter = 10 * d

  walk.list <- mclapply(1:num.iter, function (trial) {
    set.seed(seed.hash(d, theta, cut.off.mult, trial))

    raw.trace <- gen.urn.walk(d, cut.off.mult, theta)
    list(
      trial=trial,
      init.val=as.numeric(head(raw.trace, 1)),
      final.val=as.numeric(tail(raw.trace, 1))
    )
  }, mc.cores = NUM.WORKERS)

  Reduce(rbind.data.frame, walk.list)
})

write.table(urn.walks, file="./data/urn_walks.tsv", sep="\t", quote=F, row.names=F)
