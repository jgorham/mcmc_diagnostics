# Import the functions to generate an urn

if (!exists('gen.beta.binom.walk')) {

.take.step.beta.binom <- function (n, pnow) {
  # sample theta first
  # theta | x ~ Beta(x + 1, n - x + 1)
  xnow <- pnow[1]

  theta.next <- rbeta(1, xnow + 1, n - xnow + 1)
  # x | theta ~ Binom(n; theta)
  x.next <- rbinom(1, n, theta.next)
  c(x=x.next, theta=theta.next)
}

# Here l is the number of (biariate) steps to take
gen.beta.binom.walk <- function (n, l) {
  # we've gotta start from n as its hard to evaluate
  # the eigenfunctions at other spots
  pnow <- c(x=n, theta=runif(1))
  trace.walk <- as.data.frame(t(pnow))

  for (ix in 2:l) {
    pnext <- .take.step.beta.binom(n, pnow)
    trace.walk <- rbind(trace.walk, pnext)
    pnow <- pnext
  }
  trace.walk
}

}
