# -*- coding: utf-8 -*-
__doc__ = '''Fast hierarchical clustering routines for R and Python

Copyright © 2011 Daniel Müllner
<http://math.stanford.edu/~muellner>

This module provides fast hierarchical clustering routines. It is designed
to provide a replacement for the “linkage” function and its siblings in the
scipy.cluster.hierarchy module. You may use the methods in this module with
the same syntax as the corresponding SciPy functions. Please refer to the
SciPy documentation for further information on the input and output of these
functions.
'''

__all__ = ['single', 'complete', 'average', 'weighted', 'ward', 'centroid', 'median', 'linkage']

from sys import stderr
from numpy import ndarray, double, empty, require
from scipy.spatial.distance import is_valid_y, num_obs_y, pdist
from _fastcluster import linkage_wrap

def single(D):
    '''Single linkage clustering (alias). See the help on the “linkage”
function for further information.'''
    return linkage(D, method='single')

def complete(D):
    '''Complete linkage clustering (alias). See the help on the “linkage”
function for further information.'''
    return linkage(D, method='complete')

def average(D):
    '''Hierarchical clustering with the “average” distance update formula
(alias). See the help on the “linkage” function for further information.'''
    return linkage(D, method='average')

def weighted(D):
    '''Hierarchical clustering with the “weighted” distancqe update formula
(alias). See the help on the “linkage” function for further information.'''
    return linkage(D, method='weighted')

def ward(D):
    '''Hierarchical clustering with the “Ward” distance update formula (alias).
See the help on the “linkage” function for further information.'''
    return linkage(D, method='ward')

def centroid(D):
    '''Hierarchical clustering with the “centroid” distance update formula
(alias). See the help on the “linkage” function for further information.'''
    return linkage(D, method='centroid')

def median(D):
    '''Hierarchical clustering with the “median” distance update formula
(alias). See the help on the “linkage” function for further information.'''
    return linkage(D, method='median')

mthidx = {'single'   : 0,
          'complete' : 1,
          'average'  : 2,
          'weighted' : 3,
          'ward'     : 4,
          'centroid' : 5,
          'median'   : 6 }

def linkage(D, method='single', metric='euclidean', preserve_input=True):
    '''Hierarchical (agglomerative) clustering on a dissimilarity matrix or on
Euclidean data.

The argument D is either a compressed distance matrix or a collection
of m observation vectors in n dimensions as an (m×n) NumPy array. Apart
from the argument preserve_input, the methods have the same input
parameters and output format as the functions of the same name in the
package scipy.cluster.hierarchy. Therefore, the documentation is not
duplicated here. Please refer to the SciPy documentation for further
details.

The additional, optional argument preserve_input specifies whether the
fastcluster package first copies the distance matrix or writes into
the existing array. If the distance matrix is only generated for the
clustering step and is not needed afterwards, half the memory can be
saved by specifying preserve_input=False.

Note that the input array D contains unspecified values after this
procedure. It is therefore safer to write

    linkage(D, method="…", preserve_distance=False)
    del D

to make sure the matrix D is not accidentally used after it has been
used as scratch memory.

The single linkage algorithm does not write to the distance matrix or
its copy anyway, so the preserve_distance flag has no effect in this
case.'''
    if not isinstance(D, ndarray):
        raise ValueError('The first argument must be of type numpy.ndarray.')
    if len(D.shape)==1:
        if method=='single':
            assert D.dtype==double
            D_ = require(D, dtype=double, requirements=['C'])
            if D_ is not D:
                stderr.write('The condensed distance matrix had to be copied since it has the following flags:\n')
                stderr.write(str(D.flags) + '\n')
        elif preserve_input:
            D_ = D.copy()
            assert D_.dtype == double
            assert D_.flags.c_contiguous
            assert D_.flags.owndata
            assert D_.flags.writeable
            assert D_.flags.aligned
        else:
            assert D.dtype==double
            D_ = require(D, dtype=double, requirements=['C', 'A', 'W', 'O'])
            if D_ is not D:
                stderr.write('The condensed distance matrix had to be copied since it has the following flags:\n')
                stderr.write(str(D.flags) + '\n')

        is_valid_y(D_, throw=True)

        N = num_obs_y(D_)
        Z = empty((N-1,4))
        if N > 1:
            linkage_wrap(N, D_, Z, mthidx[method])
        return Z
    else:
        assert len(D.shape)==2
        N = D.shape[0]
        Z = empty((N-1,4))
        D_ = pdist(D, metric)
        assert D_.dtype == double
        assert D_.flags.c_contiguous
        assert D_.flags.owndata
        assert D_.flags.writeable
        assert D_.flags.aligned
        if N > 1:
            linkage_wrap(N, D_, Z, mthidx[method])
        return Z
