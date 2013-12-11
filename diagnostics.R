# Them diagnostics

if (!exists('.CUSUMDiag')) {

library(plyr)
library(coda)

# This diagnostic intends to find the burn-in period.
# However this is only for a quantile of some statistic.
# So instead we'll try the 10th, 20th, ... 90th percentiles
# and pick the largest burn in time.
.rafteryLewisDiag <- function (whole.trace, Qs=seq(0.1, 0.9, by=0.1)) {
  r <- 0.02    # precision of quantile estimate
  s <- 0.95    # prob we land within r of the quantile
  num.param <- ncol(whole.trace)

  quant.res <- adply(Qs, 1, function (q) {
    rl.diag <- raftery.diag(whole.trace, q=q, r=r, s=s)

    if (is.character(rl.diag$resmatrix)) {
      mn.res <- data.frame(
        M=rep(NA, num.param),
        N=rep(NA, num.param)
      )
      param.names <- colnames(whole.trace)
    } else {
      param.names <- rownames(rl.diag$resmatrix)
      mn.res <- as.data.frame(rl.diag$resmatrix[, c(1,2)])
    }

    cbind(
      param=param.names,
      q=q,
      mn.res
    )
  })[, -1]

  ddply(quant.res, .(param), function (df) {
    if (sum(is.na(df$M)) > 0) {
      max.burnin <- NA
      converged <- NA
    } else {
      max.burnin <- max(df$M)
      converged <- (max.burnin < nrow(whole.trace))
    }
    c(
      heuristic="RafteryLewis",
      stat=max.burnin,
      stat.type="max.burnin",
      has.converged=converged
    )
  })
}

.gelmanDiag <- function (trace.list) {
  g.diag <- gelman.diag(trace.list)

  data.frame(
    param=rownames(g.diag$psrf),
    heuristic="Gelman",
    stat=as.numeric(g.diag$psrf[, 1]),
    stat.type="PSRF",
    has.converged=(as.numeric(g.diag$psrf[, 1]) <= 1.1)
  )
}

.gewekeDiag <- function (half.trace) {
  g.diag <- geweke.diag(half.trace)
  geweke.stat <- g.diag$z
  geweke.conv <- (2 * pnorm(abs(geweke.stat), lower.tail=F) >= 0.05)

  data.frame(
    param=names(geweke.stat),
    heuristic="Geweke",
    stat=as.numeric(geweke.stat),
    stat.type="z.score",
    has.converged=as.logical(geweke.conv)
  )
}

.heidelDiag <- function (half.trace) {
  h.diag <- heidel.diag(half.trace)

  heidel.stat <- h.diag[, 3]  # pvalue
  heidel.conv <- h.diag[, 1] == 1
  data.frame(
    param=rownames(h.diag),
    heuristic="Heidel",
    stat=heidel.stat,
    stat.type="p.val",
    has.converged=heidel.conv
  )
}

.CUSUMDiag <- function(half.raw.trace) {
  res <- adply(half.raw.trace, 2, function (param.trace) {
    param.trace <- as.numeric(param.trace[,1])
    mu.hat <- mean(param.trace)
    signs <- sign(param.trace - mu.hat)
    D <- sum(abs(diff(signs))) / 2
    # D should be Binomial prob=1/2 and N - 1 flips
    q <- pbinom(D, length(param.trace) - 1, prob=0.5)
    pval <- (1 - 2 * abs(q - 0.5))

    data.frame(
      heuristic="CUSUM",
      stat=pval,
      stat.type="p.val",
      has.converged=(pval > 0.05)
    )
  })[, -1]

  cbind(
    param=colnames(half.raw.trace),
    res
  )
}

}
