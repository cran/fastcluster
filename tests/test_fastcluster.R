#  fastcluster: Fast hierarchical clustering routines for R and Python
#
#  Copyright © 2011 Daniel Müllner
#  <http://math.stanford.edu/~muellner>

seed = as.integer(runif(1, 0, 1e9))
set.seed(seed)

print_seed <- function() {
  return(sprintf('
Please send a bug report to the author of the \'fastcluster\' package,
Daniel Müllner. For contact details, see <http://math.stanford.edu/~muellner>.
To make the error reproducible, you must include the following number (the
random seed value) in your error report: %d.\n\n', seed))
}

compare <- function(dg1, dg2) {
  # If there is a tie (two merging steps at the same height),
  # skip this test - the situation is just too complicated
  # for this simple test file.
  h1 <- dg1$height
  h2 <- dg2$height
  rdiff1 = diff(h1) / h1[2:length(h2)]
  rdiff2 = diff(h2) / h2[2:length(h2)]
  mindiff1 = min(abs(rdiff2))
  mindiff2 = min(abs(rdiff1))
  if (max(mindiff1,mindiff2)<1e-14) {
    cat('Tie: two nearly equal merging heights. Skip this case.\n')
    return (TRUE)
  }
  # "merge" matrices must be identical
  if (!identical(dg1$merge,dg2$merge)) {
    cat('Merge matrices differ!\n')
    return(FALSE)
  }
  # "height" vectors may have small numerical errors.
  rel_error <-  max(abs(h1-h2)/pmax(abs(h1),abs(h2)))
  # We allow a relative error of 1e-14.
  if (rel_error>1e-14) {
    cat(sprintf('Height vectors differ! The maximum relative error is %e.\n', rel_error))
    return(FALSE)
  }
  if (!identical(dg1$order,dg2$order)) {
    cat('Order vectors differ!\n')
    return(FALSE)
  }
  return(TRUE)
}


test_uniform <- function() {
  n = sample(10:1000,1)
  range_exp = runif(1,min=-10, max=10)
  cat(sprintf("Number of sample points: %d\n",n))
  cat(sprintf("Dissimilarity range: [0,%g]\n",10^range_exp))
  d = runif(n*(n-1)/2, min=0, max=10^range_exp)
  # Fake a compressed distance matrix
  attributes(d) <- NULL
  attr(d,"Size") <- n
  attr(d, "call") <- 'N/A'
  class(d) <- "dist"
  test(d)
}

test_eucl <- function() {
  n = sample(10:1000,1)
  dim = sample(2:20,1)

  cat (sprintf("Number of sample points: %d\n",n))
  cat (sprintf("Dimension: %d\n",dim))

  pcd = matrix(rnorm(n*dim), c(n,dim))
  d = dist(pcd)
  test(d)
}

test <-  function(d) {
  for (method in c('single','complete','average','mcquitty','ward','centroid','median') ) {
    cat(paste('Method :', method, '\n'))
    dg_stats      = stats::hclust(d, method=method)
    dg_fastcluster = fastcluster::hclust(d, method=method)
    if (!compare(dg_stats, dg_fastcluster)) {
      stop(print_seed())
    }
  }
  cat('Passed.\n')
}

N = 25
for (i in (1:N)) {
  if (i%%2==1) {
    cat(sprintf('Random test %d of %d (uniform distribution of distances):\n',i,N))
    test_uniform()
  }
  else {
    cat(sprintf('Random test %d of %d (Gaussian density):\n',i,N))
    test_eucl()
  }
}
cat('Done.\n')
