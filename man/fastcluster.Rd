\name{fastcluster}
\alias{fastcluster}
\alias{fastcluster-package}
\alias{hclust}
\docType{package}
\title{Fast hierarchical clustering routines for R and Python}
\description{
  This function implements hierarchical clustering with the same interface as \code{\link{hclust}} from the \code{\link{stats}} package but with much faster algorithms.
}
\usage{hclust(d, method = "complete", members=NULL)}
\arguments{
\item{d}{a dissimilarity structure as produced by \code{dist}.}
\item{method}{the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of \code{"ward"}, \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"median"} or \code{"centroid"}.}
\item{members}{\code{NULL} or a vector with length size of \code{d}.}
}
\details{When the package is loaded, it overwrites the function \code{\link{hclust}} with the new code. Please refer to the original documentation on \code{\link[stats]{hclust}} for a detailed description, since the interface and the output format are exactly the same.

See the author's home page for a performance comparison and for further documentation: \url{http://math.stanford.edu/~muellner}}
\references{\url{http://math.stanford.edu/~muellner}}
\author{Daniel Müllner}
\seealso{\code{\link{hclust}}}
\examples{# Taken from stats::hclust
#
# hclust(...)        # new method
# stats::hclust(...) # old method

require(fastcluster)
require(graphics)

hc <- hclust(dist(USArrests), "ave")
plot(hc)
plot(hc, hang = -1)

## Do the same with centroid clustering and squared Euclidean distance,
## cut the tree into ten clusters and reconstruct the upper part of the
## tree from the cluster centers.
hc <- hclust(dist(USArrests)^2, "cen")
memb <- cutree(hc, k = 10)
cent <- NULL
for(k in 1:10){
  cent <- rbind(cent, colMeans(USArrests[memb == k, , drop = FALSE]))
}
hc1 <- hclust(dist(cent)^2, method = "cen", members = table(memb))
opar <- par(mfrow = c(1, 2))
plot(hc,  labels = FALSE, hang = -1, main = "Original Tree")
plot(hc1, labels = FALSE, hang = -1, main = "Re-start from 10 clusters")
par(opar)
}
\keyword{multivariate}
\keyword{cluster}
