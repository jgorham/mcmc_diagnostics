source('seed_hash.R')
source('mcmc_urn.R')
source('diagnostics.R')

library(parallel)

NMONTE <- 1:100
GELMAN.CHAINS <- 10

NUM.WORKERS <- 8
THETAS <- c(0.25, 0.5, 1)
DIM.VALUES <- c(100, 250, 1000, 4000)
CUT.OFF.MULT <- seq(0.4, 4, by=0.2)
SEED.BASE <- '1045ww'

params <- expand.grid(
  d=DIM.VALUES,
  cut.off.mult=CUT.OFF.MULT,
  theta=THETAS
)

urn.diag <- adply(params, 1, function (param) {
  d = param$d
  cut.off.mult = param$cut.off.mult
  theta = param$theta
  print(param)

  non.gelman <- mclapply(NMONTE, function (trial) {
    set.seed(seed.hash(d, theta, cut.off.mult, trial))

    raw.trace <- gen.urn.walk(d, cut.off.mult, theta)
    ll <- nrow(raw.trace)
    half.raw.trace <- as.data.frame(raw.trace[floor(ll / 2):ll, ])
    colnames(half.raw.trace) <- colnames(raw.trace)

    whole.trace <- as.mcmc(raw.trace)
    half.trace <- as.mcmc(half.raw.trace)

    res <- rbind(
      .rafteryLewisDiag(whole.trace),
      .gewekeDiag(half.trace),
      .heidelDiag(half.trace),
      .CUSUMDiag(half.raw.trace)
    )

    res <- res[with(res, order(heuristic, param)), ]
    param.names <- as.character(sort(unique(res$param)))
    init.val <- as.numeric(raw.trace[1, param.names])
    end.val <- as.numeric(raw.trace[ll, param.names])

    cbind(
      walk.length=ll,
      trial=trial,
      init.val=init.val,
      final.val=end.val,
      res
    )
  }, mc.cores = NUM.WORKERS)

  gelman <- mclapply(NMONTE, function (trial) {
    chains <- list()

    for (i in 1:GELMAN.CHAINS) {
      set.seed(seed.hash(d, theta, cut.off.mult, trial, i))

      raw.trace <- gen.urn.walk(d, cut.off.mult, theta)
      trace.walk <- as.mcmc(raw.trace)
      ll <- nrow(raw.trace)

      chains[[i]] <- trace.walk
    }

    res <- .gelmanDiag(chains)
    cbind(
      walk.length=ll,
      trial=trial,
      init.val=NA,  # we dont need this for gelman
      final.val=NA,
      res
    )
  }, mc.cores = NUM.WORKERS)

  all.df <- rbind(
    Reduce(rbind, non.gelman),
    Reduce(rbind, gelman)
  )

  cbind(
    "dim"=d,
    all.df
  )
})

write.table(urn.diag, file="./data/urn_diag.tsv", sep="\t", quote=F, row.names=F)
