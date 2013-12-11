library(ggplot2)
library(plyr)

#########################
# Final Location Figure #
#########################
urn.walks <- read.csv(file="./data/urn_walks.tsv", sep="\t")

dim.vals <- unique(urn.walks$d)
stationary <- adply(dim.vals, 1, function (d) {
  xval <- seq(0, d, by=1);
  ps <- sapply(xval, function (xx) {
    dbinom(xx, d, prob=0.5)
  })

  data.frame(
    theta=1,
    "d"=d,
    "xval"=xval,
    "stationary.probs"=ps
  )
})[, -1]

png(file="./figs/urn_stationary_distro_600x900.png", width=600, height=900)
ggplot(data=urn.walks, aes(x=final.val)) +
  geom_histogram(aes(y=..density..), color="lightblue", fill="lightblue", binwidth=1) +
  geom_line(aes(x=xval, y=stationary.probs), data=stationary, color="black", size=I(0.75)) +
  facet_grid(cut.off.mult ~ d, scales="free") +
  labs(x="Final Location",
       y="Empirical Probability",
       title="Illustration of Urn Walk Endpoint with Uniform Prior")
dev.off()

#####################
# Example traceplot #
#####################
source('seed_hash.R')
source('mcmc_urn.R')

SEED.BASE <- '1045ww'
d <- 1000
cut.off.mult <- 2
theta <- 0.25
trial <- 3
set.seed(seed.hash(d, theta, cut.off.mult, trial))

raw.trace <- gen.urn.walk(d, cut.off.mult, theta)
ll <- nrow(raw.trace)
half.raw.trace <- as.data.frame(raw.trace[floor(ll / 2):ll, ])
colnames(half.raw.trace) <- colnames(raw.trace)

half.raw.trace <- cbind(half.raw.trace, iter=(floor(ll / 2):ll))

png(file="./figs/urn_half_traceplot_d1000_l2_th1_t3_600x900.png", width=600, height=600)
ggplot(data=half.raw.trace, aes(x=iter, y=x)) +
  geom_line(color="black") +
  geom_hline(yintercept=mean(half.raw.trace$x), color="red") +
  labs(x="Iteration Number", y="X value",
       title="(Half) Traceplot for Ehrenfest Urn")
dev.off()

#####################
# Diagnostic Figure #
#####################
urn.diag <- read.csv(file="./data/urn_diag.tsv", sep="\t")

urn.diag$d <- ordered(urn.diag$d, levels=sort(unique(urn.diag$d)))
urn.diag$heuristic <- as.factor(urn.diag$heuristic)
urn.diag$theta.lbl <- as.factor(
  sapply(urn.diag$theta, function (t) {paste('theta=', t, sep='')})
)

urn.prob.conv <- ddply(urn.diag, .(heuristic, theta.lbl, d, cut.off.mult), function (df) {
  num.non.NA <- sum(!is.na(df$has.converged))
  num.converged <- sum(!is.na(df$has.converged) & df$has.converged)

  if (num.non.NA > 0) {
    bern.test <- prop.test(num.converged, num.non.NA)
    c(
      "prob.converged"=as.numeric(bern.test$estimate),
      "prob.lo"=bern.test$conf.int[1],
      "prob.hi"=bern.test$conf.int[2],
      "prob.NA"=sum(is.na(df$has.converged))/nrow(df)
    )
  } else {
    c(
      "prob.converged"=NA,
      "prob.lo"=NA,
      "prob.hi"=NA,
      "prob.NA"=1
    )
  }
})

urn.prob.conv$cut.off.mult <- ordered(urn.prob.conv$cut.off.mult,
  levels=sort(unique(urn.prob.conv$cut.off.mult)))

png(file="./figs/urn_diagnostic_by_heuristic_600x800.png", width=600, height=800)
ggplot(data=urn.prob.conv, aes(x=cut.off.mult, y=prob.converged)) +
  geom_boxplot(aes(color=heuristic, ymin=prob.lo, ymax=prob.hi)) +
  geom_point(aes(color=heuristic), size=I(0.5)) +
  facet_grid(d ~ theta.lbl) +
  labs(x="Ratio of Chain Length to Asymptotic Cut Off",
       y="Probability of Concluding Convergence",
       title="Diagnostic Probabilities For Ehrenfest Urn") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

png(file="./figs/urn_gelman_rubin_800x800.png", width=800, height=800)
ggplot(data=subset(urn.diag, heuristic == "Gelman"), aes(x=cut.off.mult, group=cut.off.mult)) +
  geom_boxplot(aes(y=stat, color=theta.lbl)) +
  facet_grid(d ~ theta.lbl, scales="free_y") +
  labs(x="Ratio of Chain Length to Asymptotic Cut Off",
       y="PSRF",
       title="Gelman-Rubin For Ehrenfest Urn")
dev.off()

# Now for the beta-binomial
bb.diag <- read.csv(file="./data/beta_binomial_diag.tsv", sep="\t")

bb.diag$N <- ordered(bb.diag$N, levels=sort(unique(bb.diag$N)))
bb.diag$heuristic <- as.factor(bb.diag$heuristic)
bb.diag$param <- as.factor(bb.diag$param)

bb.prob.conv <- ddply(bb.diag, .(heuristic, param, N, cut.off.mult), function (df) {
  num.non.NA <- sum(!is.na(df$has.converged))
  num.converged <- sum(!is.na(df$has.converged) & df$has.converged)

  if (num.non.NA > 0) {
    bern.test <- prop.test(num.converged, num.non.NA)
    c(
      "prob.converged"=as.numeric(bern.test$estimate),
      "prob.lo"=bern.test$conf.int[1],
      "prob.hi"=bern.test$conf.int[2],
      "prob.NA"=sum(is.na(df$has.converged))/nrow(df)
    )
  } else {
    c(
      "prob.converged"=NA,
      "prob.lo"=NA,
      "prob.hi"=NA,
      "prob.NA"=1
    )
  }
})

bb.prob.conv$cut.off.mult <- ordered(bb.prob.conv$cut.off.mult,
  levels=sort(unique(bb.prob.conv$cut.off.mult)))

png(file="./figs/beta_binomial_diagnostic_by_heuristic_800x1200.png", width=800, height=1200)
ggplot(data=bb.prob.conv, aes(x=cut.off.mult, y=prob.converged)) +
  geom_boxplot(aes(color=heuristic, ymin=prob.lo, ymax=prob.hi)) +
  geom_point(aes(color=heuristic), size=I(0.5)) +
  facet_grid(N ~ param) +
  labs(x="Ratio of Chain Length to Asymptotic Cut Off",
       y="Probability of Concluding Convergence",
       title="Diagnostic Probabilities For Beta Binomial Gibbs Sampler") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

png(file="./figs/beta_binomial_gelman_rubin_800x800.png", width=800, height=800)
ggplot(data=subset(bb.diag, heuristic == "Gelman"), aes(x=cut.off.mult, group=cut.off.mult)) +
  geom_boxplot(aes(y=stat, color=param)) +
  facet_grid(N ~ param, scales="free_y") +
  labs(x="Ratio of Chain Length to Asymptotic Cut Off",
       y="PSRF",
       title="Gelman-Rubin For Beta Binomial")
dev.off()
