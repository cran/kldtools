ksboot  <- function(x, y, nboots=1000, alternative=c("two.sided", "less", "greater")) {
 tol <- sqrt(.Machine$double.eps)
 bbcount <- 0
 ks.f  <- ks.test(x, y, alternative=alternative)
 nx <- length(x)
 w <- c(x, y)
 obs <- length(w)
 for (bb in 1:nboots) {
  x.s  <- sample(1:obs, obs, replace=TRUE)
  xtmp <- w[x.s[1:nx]]
  ytmp <- w[x.s[(nx + 1):obs]]
  ks.s <- suppressWarnings(ks.test(xtmp, ytmp, alternative=alternative)$statistic)
  if (ks.s >= (ks.f$statistic - tol)) bbcount <- bbcount + 1
 }
 ksboot.pval  <- bbcount/nboots
 list(ksboot.pvalue=ksboot.pval, nboots=nboots)
}
