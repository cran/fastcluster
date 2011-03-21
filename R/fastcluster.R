#  fastcluster: Fast hierarchical clustering routines for R and Python
#
#  Copyright © 2011 Daniel Müllner
#  <http://math.stanford.edu/~muellner>

hclust <- function(d, method="complete",  members = NULL)
{
  # Hierarchical clustering, on raw input data.
  METHODS <- c("single", "complete", "average", "mcquitty", "ward", "centroid", "median")
  method <- pmatch(method, METHODS)
  if (is.na(method))
    stop("Invalid clustering method.")
  if (method == -1)
    stop("Ambiguous clustering method.")
  dendrogram <- c( .Call(fastcluster, attr(d, "Size"), method, d, members),
    list(
      labels = attr(d, "Labels")
      ,method = METHODS[method]
      ,call = match.call()
      ,dist.method = attr(d, "method")
    )
  )
  class(dendrogram) <- "hclust"
  return (dendrogram)
}
