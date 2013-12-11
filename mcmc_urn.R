# Import the functions to generate an urn

if (!exists('gen.urn.walk')) {

.take.step.urn <- function (theta, d, xnow) {
  bounds <- c(xnow / d, 1 - theta * (1 - (xnow / d)))
  stepsize <- sum(runif(1) >= bounds) - 1
  xnow + stepsize
}

gen.urn.walk <- function (d, m, theta) {
  xnow <- sample(0:d, 1)
  trace.walk <- c(xnow)

  k <- round(
    m * (d / (2 * (1 + theta))) * log(d * theta),
    0
  )

  for (ix in 2:k) {
    xnext <- .take.step.urn(theta, d, xnow)
    trace.walk <- c(trace.walk, xnext)
    xnow <- xnext
  }
  data.frame(x=trace.walk)
}

}
