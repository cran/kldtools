kldtools <- function(x, y, threshold=0.975) {

## checks
if(!identical(length(x), length(y))) stop("both vectors must have the same length")
if(any(length(x) < 2)) stop("length of both vectors must be more than one")
if(any(!is.finite(x)) || any(!is.finite(y))) stop("both vectors must have olnly finite values")
if(any(x < 0) || any(y < 0)) stop("both vectors must consist of nonnegative values")
if(!all.equal(y, as.integer(y)) || !all.equal(y, as.integer(y))) stop("both vectors must consist of whole numbers")

## The numerical constant: by default, the 97.5% quantile of the normal distribution
qq <- qnorm(threshold)

## The main variables
a <- x/sum(x)
b <- y/sum(y)

## The technical variables, useful as shortcuts
na <- length(a)
nb <- length(b)
na1 <- na - 1
nb1 <- nb - 1

## calculate Kullback-Leibler divergence (KLD)
KLD <- sum(a * log(a/b))

## calculate symmetrized KLD
KLD.s <- mean(c(KLD, sum(b * log(b/a))))

## handling zeros
b[b == 0] <- 1/sum(y)

## calculate the g(v) vector used for the calculation of variance
g.v1 = log((a[1:(na1)] * b[nb]) / (b[1:(nb1)] * a[na]))
g.v2 = -a[1:(na1)] / b[1:(nb1)] + (a[na] / b[nb])
g.V <- c(g.v1, g.v2)

## calculate the same vector for symmetrized KLD
g.x <- 1/2*(log(a[1:na1]/b[1:nb1]) - log(a[na]/b[nb])) - 1/2*(b[1:nb1]/a[1:na1] - b[nb]/a[na])
g.y <-  1/2*(log(b[1:nb1]/a[1:na1]) - log(b[nb]/a[na])) - 1/2*(a[1:na1]/b[1:nb1] - a[na]/b[nb])
g.S <- c(g.x, g.y)

## construct the first quasi-diagonal matrix for the calculation of variance
ma1 <- matrix(rep(a[1:(na1)], na1), nrow=na1, byrow=TRUE)
ma2 <- sweep(ma1, MARGIN=1, -a[1:(na1)], `*`)
diag(ma2) <- a[1:(na1)]*(1-a[1:(na1)])

## construct the second quasi-diagonal matrix for the calculation of variance
mb1 <- matrix(rep(b[1:(nb1)], nb1), nrow=nb1, byrow=TRUE)
mb2 <- sweep(mb1, MARGIN=1, -b[1:(nb1)], `*`)
diag(mb2) <- b[1:(nb1)]*(1-b[1:(nb1)])
lambda <- sum(x)/sum(y)
mb3 <- lambda * mb2

## combine both matrices from the above
m0 <- matrix(0, nrow=na1, ncol=na1)
mm <- rbind(cbind(ma2, m0), cbind(m0, mb3))

## calculate the variance using the g(v) vector and the combined matrix
mu <- mm %*% g.V
sd <- c(sqrt(t(g.V) %*% mu))

## calculate the variance for symmetrized KLD
mu.s <- mm %*% g.S
sd.s <- c(sqrt(t(g.S) %*% mu.s))

## calculate the lower and upper limits of the confidence interval
left <- KLD - (qq * sd/sqrt(sum(x)))
right <- KLD + (qq * sd/sqrt(sum(x)))

## calculate the lower and upper limits of the confidence interval from symmetrized KLD
left.s <- KLD.s - (qq * sd.s/sqrt(sum(x)))
right.s <- KLD.s + (qq * sd.s/sqrt(sum(x)))

## calculate Turing's perspective estimator
.Turing <- function(a, b) {
va <- sum(a) - a
vb <- sum(b) - b
sumb <- suma <- numeric(length(a))
for (k in 1:length(a)) { # k-part
 ## left half
 blong <- 1 - b[k]/(sum(b) - 1:vb[k] + 1) # j-part
 sumb[k] <- sum(1/(1:vb[k]) * cumprod(blong[1:vb[k]])) # v-part
 ## right half
 along <- 1 - (a[k] - 1)/(sum(a) - 1:va[k])
 suma[k] <- sum(1/(1:va[k]) * cumprod(along[1:va[k]]))
 }
sum(a/sum(a) * (sumb - suma))
}
Turing <- .Turing(x, y)

## calculate the lower and upper limits of the confidence interval from Turing's estimator
left.Turing <- abs(Turing) - (qq * sd/sqrt(sum(x)))
right.Turing <- abs(Turing) + (qq * sd/sqrt(sum(x)))

## return results as a list of values
list(KLD=KLD, KLD.sd=sd, KLD.ci.left=left, KLD.ci.right=right,
 KLD.s=KLD.s, KLD.s.sd=sd.s, KLD.s.ci.left=left.s, KLD.s.ci.right=right.s,
 Turing=Turing, Turing.sd=sd, Turing.left=left.Turing, Turing.right=right.Turing)
}
