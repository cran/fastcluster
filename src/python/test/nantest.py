'''Test whether the fastcluster package correctly recognizes NaN values
and raises a FloatingPointError.'''
import numpy as np
import fastcluster

n = np.random.random_integers(2,100)

# Part 1: distance matrix input

N = n*(n-1)//2
D = np.random.rand(N)
# Insert a single NaN value
pos = np.random.randint(N)
D[pos] = np.nan

for method in ['single', 'complete', 'average', 'weighted', 'ward',
               'centroid', 'median']:
    try:
        fastcluster.linkage(D, method=method)
        raise AssertionError('fastcluster did not detect a NaN value!')
    except FloatingPointError:
        pass

# Next: the original array does not contain a NaN, but a NaN occurs
# as an updated distance.
for method in ['average', 'weighted', 'ward', 'centroid', 'median']:
    try:
        fastcluster.linkage([np.inf,-np.inf,-np.inf], method=method)
        raise AssertionError('fastcluster did not detect a NaN value!')
    except FloatingPointError:
        pass

# Part 2: vector input

dim = np.random.random_integers(2,12)
X = np.random.rand(n,dim)
pos = (np.random.randint(n), np.random.randint(dim))
# Insert a single NaN coordinate
X[pos] = np.nan

for method in ['single', 'ward', 'centroid', 'median']:
    try:
        fastcluster.linkage_vector(X, method=method)
        raise AssertionError('fastcluster did not detect a NaN value!')
    except FloatingPointError:
        pass

print('OK.')
