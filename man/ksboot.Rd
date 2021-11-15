\name{ksboot}
\alias{ksboot}

\title{Bootstrapping based on the Kolmogorov-Smirnov test}

\description{Performs bootstrapping with the Kolmogorov-Smirnov test to
estimate differences between frequencies}

\usage{

ksboot(x, y, nboots=1000, alternative=c("two.sided", "less", "greater"))

}

\arguments{

\item{x}{The first vector}

\item{y}{The second vector}

\item{nboots}{The number of bootstraps to perform}

\item{alternative}{The type of alternative hypothesis (the default is
"two.sided")}

}

\details{

This bootstrap version of the Kolmogorov-Smirnov test is suitable for
estimating not only continuous but also frequency distributions. This is
because bootstrap theories suggest that the asymptotic theory of
estimates (which is built on the bootstrapping data) is, in a sense,
similar to the asymptotic theory of large data sets. Hence, although the
Kolmogorov-Smirnov test is initially designed for continuous
distributions, in bootstrapping, it is possible to apply this method to
discrete random variables, for which an empirical distribution function
is built on the observed frequencies (see Abadie 2002 for an example).

}

\value{

The list with the following components: "ksboot.pvalue" for the bootstrap
p-value of the Kolmogorov-Smirnov test, calculated for the null hypothesis
that the probability densities of two compared distributions are the
same; "nboots" for the number of the completed bootstraps.

}

\references{

Abadie A. 2002. Bootstrap tests for distributional treatment effects in
instrumental variable models. Journal of the American statistical
Association. 97: 284-292.

}

\author{Jasjeet S. Sekhon, Alexey Shipunov}

\seealso{\code{\link{ks.test}}}

\examples{

data <- stack(data.frame(V1=c(1213, 57683, 74466, 44419, 17481, 3403, 42252, 7045,
 29445, 15004, 21337, 1892, 21861, 238, 26574, 17579),
 V2=c(3185, 29692, 12570, 26081, 4992, 1659, 16592, 1748, 37583, 6751, 10188, 355,
 8116, 9, 5064, 1846)))

ksboot(data$values[data$ind == "V1"], data$values[data$ind == "V2"])

pairwise.table(function(i, j)
 suppressWarnings(ksboot(data$values[as.integer(data$ind) == i],
 data$values[as.integer(data$ind) == j])$ksboot.pvalue),
 levels(data$ind), p.adjust.method="bonferroni")

pairwise.table(function(i, j)
 suppressWarnings(ksboot(data$values[as.integer(data$ind) == i],
 data$values[as.integer(data$ind) == j])$ksboot.pvalue),
 levels(data$ind), p.adjust.method="none")

}

\keyword{htest}
