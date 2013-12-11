source('seed_hash.R')
source('mcmc_betabinom.R')
source('diagnostics.R')

library(parallel)

NMONTE <- 1:100
GELMAN.CHAINS <- 10

NUM.WORKERS <- 3
N <- c(100, 250, 1000, 2000, 4000)
CUT.OFF.MULT <- seq(0.4, 4, by=0.2)
SEED.BASE <- 'ljuh234'

params <- expand.grid(
  N=N,
  cut.off.mult=CUT.OFF.MULT
)

bb.diag <- adply(params, 1, function (param) {
  N = param$N
  cut.off.mult = param$cut.off.mult
  k <- round(N * cut.off.mult)
  print(param)

  non.gelman <- mclapply(NMONTE, function (trial) {
    set.seed(seed.hash(N, cut.off.mult, trial))

    raw.trace <- gen.beta.binom.walk(N, k)
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
      set.seed(seed.hash(N, cut.off.mult, trial, i))

      raw.trace <- gen.beta.binom.walk(N, k)
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
    "N"=N,
    all.df
  )
})

write.table(bb.diag, file="./data/beta_binomial_diag.tsv", sep="\t", quote=F, row.names=F)
