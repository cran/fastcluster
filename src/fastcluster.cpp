/*
  fastcluster: Fast hierarchical clustering routines for R and Python

  Copyright © 2011 Daniel Müllner
  <http://math.stanford.edu/~muellner>
*/
#include <float.h>
#ifndef DBL_MANT_DIG
#error The constant DBL_MANT_DIG could not be defined.
#endif

#include <math.h>
#include <algorithm>

#ifndef LONG_MAX
#include <limits.h>
#endif
#ifndef LONG_MAX
#error The constant LONG_MAX could not be defined.
#endif
#ifndef INT_MAX
#error The constant INT_MAX could not be defined.
#endif

#ifndef INT32_MAX
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#endif
#ifndef INT32_MAX
#error The constant INT32_MAX could not be defined.
#endif

typedef int_fast32_t t_index;
#define MAX_INDEX INT32_MAX
#if (LONG_MAX < MAX_INDEX)
#error The integer format "t_index" must not have a greater range than "long int".
#endif
#if (INT_MAX > MAX_INDEX)
#error The integer format "int" must not have a greater range than "t_index".
#endif
typedef int_fast64_t long_index;
typedef double t_float;
#define T_FLOAT_MANT_DIG DBL_MANT_DIG

enum {
  // non-Euclidean methods
  METHOD_METR_SINGLE           = 0,
  METHOD_METR_COMPLETE         = 1,
  METHOD_METR_AVERAGE          = 2,
  METHOD_METR_WEIGHTED         = 3,
  METHOD_METR_WARD             = 4,
  METHOD_METR_CENTROID         = 5,
  METHOD_METR_MEDIAN           = 6
};

enum {
   // Return values
  RET_SUCCESS        = 0,
  RET_MEMORY_ERROR   = 1,
  RET_STL_ERROR      = 2,
  RET_UNKNOWN_ERROR  = 3
 };

// self-destructing array pointer
template <typename type>
class auto_array_ptr{
private:
  type * ptr;
public:
  auto_array_ptr() { ptr = NULL; }
  template <typename index>
  auto_array_ptr(index size) { init(size); }
  template <typename index, typename value>
  auto_array_ptr(index size, value val) { init(size, val); }
  ~auto_array_ptr() { delete [] ptr; }
  void free() {
    delete [] ptr;
    ptr = NULL;
  }
  template <typename index>
  void init(index size) { ptr = new type [size]; }
  template <typename index, typename value>
  void init(index size, value val) {
    init(size);
    for (index i=0; i<size; i++) ptr[i] = val;
  }
  operator type *() { return ptr; }
};

struct node {
  t_index node1, node2;
  t_float dist;

  /*
  inline bool operator< (const node a) const {
    return this->dist < a.dist;
  }
  */

  inline friend bool operator< (const node a, const node b) {
    return a.dist < b.dist;
  }
};

class cluster_result {
private:
  auto_array_ptr<node> Z;
  t_index pos;

public:
  cluster_result(const t_index size) {
    Z.init(size);
    pos = 0;
  }

  void append(t_index node1, t_index node2, t_float dist) {
    Z[pos].node1 = node1;
    Z[pos].node2 = node2;
    Z[pos].dist  = dist;
    pos++;
  }

  node * operator[] (t_index idx) { return Z + idx; }

  void sqrt() {
    for (t_index i=0; i<pos; i++) {
      Z[i].dist = ::sqrt(Z[i].dist);
    }
  }
};

class doubly_linked_list {
  /*
    Class for a doubly linked list. Initially, the list is the integer range
    [0, size]. We provide a forward iterator and a method to delete an index
    from the list.

    Typical use: for (i=L.start; L<size; i=L.succ[I])
    or
    for (i=somevalue; L<size; i=L.succ[I])
  */
public:
  t_index start;
  auto_array_ptr<t_index> succ;

private:
  auto_array_ptr<t_index> pred;
  // Not necessarily private, we just do not need it in this instance.

public:
  doubly_linked_list(const t_index size) {
    // Initialize to the given size.
    pred.init(size+1);
    succ.init(size+1);
    for (t_index i=1; i<size+1; i++)
      pred[i] = i-1;
    // pred[0] is never accessed!
    for (t_index i=0; i<size; i++)
      succ[i] = i+1;
    //succ[size] is never accessed!
    start = 0;
  }

  void remove(t_index idx) {
    // Remove an index from the list.
    if (idx==start) {
      start = succ[idx];
    }
    else {
      succ[pred[idx]] = succ[idx];
      pred[succ[idx]] = pred[idx];
    }
  }
};

// Indexing functions
// D is the upper triangular part of a symmetric  N x N matrix
// We require r_ > c_ !
#define D_(r_,c_) ( D[(static_cast<long_index>(2*N-3-(r_))*(r_)>>1)-1+(c_)] )
// Z is a (N-1) x 4 array
#define Z_(_r, _c) (Z[(_r)*4 + (_c)])

/*
 * Lookup function for a union-find data structure.
 *
 * The function finds the root of idx by going iteratively through all
 * parent elements until a root is found. An element i is a root if
 * nodes[i] is zero. To make subsequent searches faster, the entry for
 * idx and all its parents is updated with the root element.
 */
class union_find {
private:
  auto_array_ptr<t_index> parent;
  t_index nextparent;

public:
  void init(t_index size) {
    parent.init(2*size-1, 0);
    nextparent = size;
  }

  t_index Find(t_index idx) {
    if (parent[idx] !=0 ) { // a → b
      t_index p = idx;
      idx = parent[idx];
      if (parent[idx] !=0 ) { // a → b → c
        do {
          idx = parent[idx];
        } while (parent[idx] != 0);
        do {
          t_index tmp = parent[p];
          parent[p] = idx;
          p = tmp;
        } while (parent[p] != idx);
      }
    }
    return idx;
  }

  void Union (t_index node1, t_index node2) {
    parent[node1] = parent[node2] = nextparent++;
  }
};

static void single_linkage(const t_index N, const t_float * const D, cluster_result & Z2) {
/*
    N: integer
    D: condensed distance matrix N*(N-1)/2
    Z2: output data structure

    The basis of this algorithm is Rohlf's Algorithm 76:

    F.James Rohlf, “Algorithm 76,” The Computer Journal, vol. 16, 1973,
    p. 93–95.
*/
  t_index i;
  t_index idx2;
  doubly_linked_list AR(N);
  auto_array_ptr<t_float> d(N);

  t_index new_node;
  t_index prev_node;
  t_float min;

  // first iteration
  idx2 = 1;
  min = d[1] = D[0];
  for (i=2; i<N; i++) {
    d[i] = D[i-1];
    if (d[i] < min) {
      min = d[i];
      idx2 = i;
    }
  }
  Z2.append(0, idx2, min);

  for (new_node=1; new_node<N-1; new_node++) {
    prev_node = idx2;
    AR.remove(prev_node);

    idx2 = AR.succ[0];
    min = d[idx2];
    for (i=idx2; i<prev_node; i=AR.succ[i]) {
      if (d[i] > D_(i, prev_node))
        d[i] = D_(i, prev_node);
      if (d[i] < min) {
        min = d[i];
        idx2 = i;
      }
    }
    for (; i<N; i=AR.succ[i]) {
      if (d[i] > D_(prev_node, i))
        d[i] = D_(prev_node, i);
      if (d[i] < min) {
        min = d[i];
        idx2 = i;
      }
    }
    Z2.append(prev_node, idx2, min);
  }
}


inline static void f_single( t_float * b, t_float a ) {
  if (*b > a) *b = a;
}
inline static void f_complete( t_float * b, t_float a ) {
  if (*b < a) *b = a;
}
inline static void f_average( t_float * b, t_float a, t_float s, t_float t) {
  *b = s*a + t*(*b);
}
inline static void f_weighted( t_float * b, t_float a) {
  *b = (a+*b)/2;
}
inline static void f_ward( t_float * b, t_float a, t_float c, t_float s, t_float t, t_float v) {
  *b = ( (v+s)*a - v*c + (v+t)*(*b) ) / (s+t+v);
  //*b = a+(*b)-(t*a+s*(*b)+v*c)/(s+t+v);
}
inline static void f_centroid( t_float * b, t_float a, t_float c, t_float s, t_float t) {
  *b = s*a + t*(*b) - s*t*c;
}
inline static void f_median( t_float * b, t_float a, t_float c) {
  *b = (a+(*b))/2 - c/4;
}


template <const unsigned char method, typename t_members>
static void NN_chain_linkage(const t_index N, t_float * const D, t_members * const members, cluster_result & Z2) {
/*
    N: integer
    D: condensed distance matrix N*(N-1)/2
    Z2: output data structure

    This is the NN-chain algorithm, described on page 86 in the following book:

﻿   Fionn Murtagh, Multidimensional Clustering Algorithms,
    Vienna, Würzburg: Physica-Verlag, 1985.
*/
  t_index i;

  auto_array_ptr<t_index> NN_chain(N);
  t_index NN_chain_tip = 0;

  t_index idx1, idx2;

  t_float size1, size2;
  doubly_linked_list AR(N);

  t_index new_node;
  t_float min;

  for (new_node=0; new_node<N-1; new_node++) {
    if (NN_chain_tip <= 3) {
      NN_chain[0] = idx1 = AR.start;
      NN_chain_tip = 1;

      idx2 = AR.succ[idx1];
      min = D_(idx1,idx2);
      for (i=AR.succ[idx2]; i<N; i=AR.succ[i]) {
        if (D_(idx1,i) < min) {
          min = D_(idx1,i);
          idx2 = i;
        }
      }
    }  // a: idx1   b: idx2
    else {
      NN_chain_tip -= 3;
      idx1 = NN_chain[NN_chain_tip-1];
      idx2 = NN_chain[NN_chain_tip];
      min = idx1<idx2 ? D_(idx1,idx2) : D_(idx2,idx1);
    }  // a: idx1   b: idx2

    do {
      NN_chain[NN_chain_tip] = idx2;

      for (i=AR.start; i<idx2; i=AR.succ[i]) {
        if (D_(i,idx2) < min) {
          min = D_(i,idx2);
          idx1 = i;
        }
      }
      for (i=AR.succ[idx2]; i<N; i=AR.succ[i]) {
        if (D_(idx2,i) < min) {
          min = D_(idx2,i);
          idx1 = i;
        }
      }

      idx2 = idx1;
      idx1 = NN_chain[NN_chain_tip++];

    } while (idx2 != NN_chain[NN_chain_tip-2]);

    Z2.append(idx1, idx2, min);

    if (idx1>idx2) {
      t_index tmp = idx1;
      idx1 = idx2;
      idx2 = tmp;
    }

    if (method==METHOD_METR_AVERAGE ||
        method==METHOD_METR_WARD) {
      size1 = static_cast<t_float>(members[idx1]);
      size2 = static_cast<t_float>(members[idx2]);
      members[idx2] += members[idx1];
    }

    // Remove the smaller index from the valid indices (AR).
    AR.remove(idx1);

    switch (method) {
      case METHOD_METR_SINGLE: {
        /*
        Single linkage.

        Characteristic: new distances are never longer than the old distances.
        */
        // Update the distance matrix in the range [start, idx1).
        for (i=AR.start; i<idx1; i=AR.succ[i])
          f_single(&D_(i, idx2), D_(i, idx1) );
        // Update the distance matrix in the range (idx1, idx2).
        for (; i<idx2; i=AR.succ[i])
          f_single(&D_(i, idx2), D_(idx1, i) );
        // Update the distance matrix in the range (idx2, N).
        for (i=AR.succ[idx2]; i<N; i=AR.succ[i])
          f_single(&D_(idx2, i), D_(idx1, i) );
        break;
      }

      case METHOD_METR_COMPLETE: {
        /*
        Complete linkage.

        Characteristic: new distances are never shorter than the old distances.
        */
        // Update the distance matrix in the range [start, idx1).
        for (i=AR.start; i<idx1; i=AR.succ[i])
          f_complete(&D_(i, idx2), D_(i, idx1) );
        // Update the distance matrix in the range (idx1, idx2).
        for (; i<idx2; i=AR.succ[i])
          f_complete(&D_(i, idx2), D_(idx1, i) );
        // Update the distance matrix in the range (idx2, N).
        for (i=AR.succ[idx2]; i<N; i=AR.succ[i])
          f_complete(&D_(idx2, i), D_(idx1, i) );
        break;
      }

      case METHOD_METR_AVERAGE: {
        /*
        Average linkage.

        Shorter and longer distances can occur.
        */
        // Update the distance matrix in the range [start, idx1).
        t_float s = size1/(size1+size2);
        t_float t = size2/(size1+size2);
        for (i=AR.start; i<idx1; i=AR.succ[i])
          f_average(&D_(i, idx2), D_(i, idx1), s, t );
        // Update the distance matrix in the range (idx1, idx2).
        for (; i<idx2; i=AR.succ[i])
          f_average(&D_(i, idx2), D_(idx1, i), s, t );
        // Update the distance matrix in the range (idx2, N).
        for (i=AR.succ[idx2]; i<N; i=AR.succ[i])
          f_average(&D_(idx2, i), D_(idx1, i), s, t );
        break;
      }

      case METHOD_METR_WEIGHTED: {
        /*
        Weighted linkage.

        Shorter and longer distances can occur.
        */
        // Update the distance matrix in the range [start, idx1).
        for (i=AR.start; i<idx1; i=AR.succ[i])
          f_weighted(&D_(i, idx2), D_(i, idx1) );
        // Update the distance matrix in the range (idx1, idx2).
        for (; i<idx2; i=AR.succ[i])
          f_weighted(&D_(i, idx2), D_(idx1, i) );
        // Update the distance matrix in the range (idx2, N).
        for (i=AR.succ[idx2]; i<N; i=AR.succ[i])
          f_weighted(&D_(idx2, i), D_(idx1, i) );
        break;
      }

      case METHOD_METR_WARD: {
        /*
        Ward linkage.

        Shorter and longer distances can occur, not smaller than min(d1,d2)
        but maybe bigger than max(d1,d2).
        */
        // Update the distance matrix in the range [start, idx1).
        for (i=AR.start; i<idx1; i=AR.succ[i])
          f_ward(&D_(i, idx2), D_(i, idx1), min,
                 size1, size2, static_cast<t_float>(members[i]) );
        // Update the distance matrix in the range (idx1, idx2).
        for (; i<idx2; i=AR.succ[i])
          f_ward(&D_(i, idx2), D_(idx1, i), min,
                 size1, size2, static_cast<t_float>(members[i]) );
        // Update the distance matrix in the range (idx2, N).
        for (i=AR.succ[idx2]; i<N; i=AR.succ[i])
          f_ward(&D_(idx2, i), D_(idx1, i), min,
                 size1, size2, static_cast<t_float>(members[i]) );
        break;
      }
    }
  }
}

class binary_min_heap {
  /*
  Class for a binary min-heap. The data resides in an array A. The elements of A
  are not changed but two lists I and R of indices are generated which point to
  elements of A and backwards.

  The heap tree structure is

     H[2*i+1]     H[2*i+2]
         \            /
          \          /
           ≤        ≤
            \      /
             \    /
              H[i]

  where the children must be less or equal than their parent. Thus, H[0] contains
  the minimum. The lists I and R are made such that H[i] = A[I[i]] and R[I[i]] = i.
  */
private:
  t_float * A;
  t_index size;
  auto_array_ptr<t_index> I;
  auto_array_ptr<t_index> R;
public:
  binary_min_heap(const t_index size) { // size must be strictly smaller than MAX_INDEX
    // Allocate memory and initialize the lists I and R to the identity. This does
    // not make it a heap. Call heapify afterwards!
    this->size = size;
    I.init(size);
    R.init(size);
    for (t_index i=0; i<size; i++)
      R[i] = I[i] = i;
  }

  void heapify(t_float * const A) {
    // Arrange the indices I and R so that H[i] := A[I[i]] satisfies the heap
    // condition H[i] < H[2*i+1] and H[i] < H[2*i+2] for each i.
    t_index idx1, idx2;
    this->A = A;
    for (t_index i=1; i<size; i++) {
      for (idx1=i; idx1>0; idx1=idx2) {
        idx2 = (idx1-1)>>1;
        if (H(idx1)<H(idx2))
          heap_swap(idx1, idx2);
        else break;
      }
    }
  }

  inline operator t_index() {
    // Return the minimal element.
    return I[0];
  }

  void heap_pop() {
    // Remove the minimal element from the heap.
    t_index i, j;
    for (i=0; (j=2*i+1)<size;) {
      if (is_empty(j)) {
        if ( is_invalid(++j) ) break ;
      }
      else if ( !safe_geq(j+1,j) ) j++;
      heap_swap(i, j);
      i = j;
    }
    invalidate(i);
  }

  inline void update ( const t_index idx, const t_float val ) {
    // Update the element A[i] with val and re-arrange the indices the preserve the
    // heap condition.
    if (val<=A[idx])
      update_leq(idx, val);
    else
      update_geq(idx, val);
  }

  inline void update_leq ( const t_index idx, const t_float val ) {
    // Use this when the new value is not more than the old value.
    t_index i,j;
    A[idx] = val;
    for (i = R[idx]; (i>0) && H(i)<H(j=(i-1)>>1); i=j)
      heap_swap(i,j);
  }

  inline void update_geq ( const t_index idx, const t_float val ) {
    // Use this when the new value is not less than the old value.
    t_index i, j;
    A[idx] = val;
    for (i=R[idx]; (j=2*i+1)<size; i=j) {
      if ( is_empty(j) || H(j)>=H(i) ) {
        if ( safe_geq(++j,i) ) break;
      }
      else if ( !safe_geq(j+1,j) ) j++;
      heap_swap(i, j);
    }
  }

private:
  inline void heap_swap(const t_index i, const t_index j) {
    // Swap two indices.
    t_index tmp = I[i];
    I[i] = I[j];
    I[j] = tmp;
    R[I[i]] = i;
    R[I[j]] = j;
  }

  inline bool safe_geq(const t_index i, const t_index j) {
    return ( is_invalid(i) || H(i)>=H(j) );
  }

  inline bool is_invalid(const t_index i) {
    return ( i>=size || is_empty(i) );
  }

  inline t_float H(const t_index i) {
    return A[I[i]];
  }

  inline bool is_empty(const t_index i) {
    return I[i]==MAX_INDEX;
  }

  inline void invalidate(const t_index i) {
    I[i] = MAX_INDEX;
  }
};

template <const unsigned char method, typename t_members>
static void generic_linkage(const t_index N, t_float * const D, t_members * const members, cluster_result & Z2) {
  /*
    N: integer
    D: condensed distance matrix N*(N-1)/2
    Z2: output data structure
  */

  const t_index N_1 = N-1;
  t_index i, j; // loop variables
  t_index idx1, idx2; // row and column indices

  auto_array_ptr<t_index> n_nghbr(N_1); // array of nearest neighbors
  auto_array_ptr<t_float> mindist(N_1); // distances to the nearest neighbors
  auto_array_ptr<t_index> row_repr(N);// row_repr[i]: node number that the i-th row represents
  doubly_linked_list AR(N); // active rows
  binary_min_heap minidx(N_1); // minimum heap structure for the nearest neighbors

  t_index node1, node2;     // node numbers in the output
  t_float size1, size2;     // and their cardinalities

  t_float min; // minimum and row index for nearest-neighbor search
  t_index idx;

  for (i=0; i<N; i++)
    // Build a list of row ↔ node label assignments.
    // Initially i ↦ i
    row_repr[i] = i;

  // Initialize the minimal distances:
  // Find the nearest neighbor of each point.
  // n_nghbr[i] = argmin_{j>i} D(i,j) for i in range(N-1)
  t_float * DD = D;
  for (i=0; i<N_1; i++) {
    min = *(DD++);
    idx = j = i+1;
    while (j<N_1) {
      j++;
      if (*DD<min) {
        min = *DD;
        idx = j;
      }
      DD++;
    }
    mindist[i] = min;
    n_nghbr[i] = idx;
  }
  // Put the minimal distances into a heap structure to make the repeated global
  // minimum searches fast.
  minidx.heapify(mindist);

  // Main loop: We have N-1 merging steps.
  for (i=0; i<N_1; i++) {
    /*
      Here is a special feature that allows fast bookkeeping and updates of the
      minimal distances.

      mindist[i] stores a lower bound on the minimum distance of the point i to
      all points of higher index:

          mindist[i] ≥ min_{j>i} D(i,j)

      Normally, we have equality. However, this minimum may become invalid due to
      the updates in the distance matrix. The rules are:

      1) If mindist[i] is equal to D(i, n_nghbr[i]), this is the correct minimum
         and n_nghbr[i] is a nearest neighbor.

      2) If mindist[i] is smaller than D(i, n_nghbr[i]), this might not be the
         correct minimum. The minimum needs to be recomputed.

      3) mindist[i] is never bigger than the true minimum. Hence, we never miss the
         true minimum if we take the smallest mindist entry, re-compute the value if
         necessary (thus maybe increasing it) and looking for the now smallest
         mindist entry until a valid minimal entry is found. This step is done in the
         lines below.

      The update process for D below takes care that these rules are fulfilled. This
      makes sure that the minima in the rows D(i,i+1:)of D are re-calculated when
      necessary but re-calculation is avoided whenever possible.

      The re-calculation of the minima makes the worst-case runtime of this algorithm
      cubic in N. We avoid this whenever possible, and in most cases the runtime
      appears to be quadratic.
    */
    idx1 = minidx;
    if (method != METHOD_METR_SINGLE) {
      while ( D_(idx1, n_nghbr[idx1]) != mindist[idx1] ) {
        // Recompute the minimum mindist[idx1] and n_nghbr[idx1].
        n_nghbr[idx1] = j = AR.succ[idx1]; // exists, maximally N-1
        min = D_(idx1,j);
        for (j=AR.succ[j]; j<N; j=AR.succ[j]) {
          if (D_(idx1,j)<min) {
            min = D_(idx1,j);
            n_nghbr[idx1] = j;
          }
        }
        minidx.update_geq(idx1, min); // Update the heap with the new true minimum.
        idx1 = minidx; // And search for the (possibly different) minimal entry.
      }
    }

    minidx.heap_pop(); // Remove the current minimum from the heap.
    idx2 = n_nghbr[idx1];

    // Write the newly found minimal pair of nodes to the output array.
    node1 = row_repr[idx1];
    node2 = row_repr[idx2];

    if (method==METHOD_METR_AVERAGE ||
        method==METHOD_METR_WARD ||
        method==METHOD_METR_CENTROID) {
      size1 = static_cast<t_float>(members[idx1]);
      size2 = static_cast<t_float>(members[idx2]);
      members[idx2] += members[idx1];
    }
    Z2.append(node1, node2, mindist[idx1]);

    // Remove idx1 from the list of active indices (AR).
    AR.remove(idx1);
    // Index idx2 now represents the new (merged) node with label N+i.
    row_repr[idx2] = N+i;

    // Update the distance matrix
    switch (method) {
      case METHOD_METR_SINGLE: {
        /*
        Single linkage.

        Characteristic: new distances are never longer than the old distances.
        */
        // Update the distance matrix in the range [start, idx1).
        for (j=AR.start; j<idx1; j=AR.succ[j]) {
          f_single(&D_(j, idx2), D_(j, idx1));
          if (n_nghbr[j] == idx1)
            n_nghbr[j] = idx2;
        }
        // Update the distance matrix in the range (idx1, idx2).
        for (; j<idx2; j=AR.succ[j]) {
          f_single(&D_(j, idx2), D_(idx1, j));
          // If the new value is below the old minimum in a row, update
          // the mindist and n_nghbr arrays.
          if (D_(j, idx2)<mindist[j]) {
            minidx.update_leq(j, D_(j, idx2));
            n_nghbr[j] = idx2;
          }
        }
        // Update the distance matrix in the range (idx2, N).
        // Recompute the minimum mindist[idx2] and n_nghbr[idx2].
        if (idx2<N_1) {
          min = mindist[idx2];
          for (j=AR.succ[idx2]; j<N; j=AR.succ[j]) {
            f_single(&D_(idx2, j), D_(idx1, j) );
            if (D_(idx2, j) < min) {
              n_nghbr[idx2] = j;
              min = D_(idx2, j);
            }
          }
          minidx.update_leq(idx2, min);
        }
        break;
      }

      case METHOD_METR_COMPLETE: {
        /*
        Complete linkage.

        Characteristic: new distances are never shorter than the old distances.
        */
        // Update the distance matrix in the range [start, idx1).
        for (j=AR.start; j<idx1; j=AR.succ[j]) {
          f_complete(&D_(j, idx2), D_(j, idx1) );
          if (n_nghbr[j] == idx1)
            n_nghbr[j] = idx2;
        }
        // Update the distance matrix in the range (idx1, idx2).
        for (; j<idx2; j=AR.succ[j])
          f_complete(&D_(j, idx2), D_(idx1, j) );
        // Update the distance matrix in the range (idx2, N).
        for (j=AR.succ[idx2]; j<N; j=AR.succ[j])
          f_complete(&D_(idx2, j), D_(idx1, j) );
        break;
      }

      case METHOD_METR_AVERAGE: {
        /*
        Average linkage.

        Shorter and longer distances can occur.
        */
        // Update the distance matrix in the range [start, idx1).
        t_float s = size1/(size1+size2);
        t_float t = size2/(size1+size2);
        for (j=AR.start; j<idx1; j=AR.succ[j]) {
          f_average(&D_(j, idx2), D_(j, idx1), s, t);
          if (n_nghbr[j] == idx1)
            n_nghbr[j] = idx2;
        }
        // Update the distance matrix in the range (idx1, idx2).
        for (; j<idx2; j=AR.succ[j]) {
          f_average(&D_(j, idx2), D_(idx1, j), s, t);
          if (D_(j, idx2)<mindist[j]) {
            minidx.update_leq(j, D_(j, idx2));
            n_nghbr[j] = idx2;
          }
        }
        // Update the distance matrix in the range (idx2, N).
        if (idx2<N_1) {
          n_nghbr[idx2] = j = AR.succ[idx2]; // exists, maximally N-1
          f_average(&D_(idx2, j), D_(idx1, j), size1, size2 );
          min = D_(idx2,j);
          for (j=AR.succ[j]; j<N; j=AR.succ[j]) {
            f_average(&D_(idx2, j), D_(idx1, j), s, t);
            if (D_(idx2,j)<min) {
              min = D_(idx2,j);
              n_nghbr[idx2] = j;
            }
          }
          minidx.update(idx2, min);
        }
        break;
      }

      case METHOD_METR_WEIGHTED: {
        /*
        Weighted linkage.

        Shorter and longer distances can occur.
        */
        // Update the distance matrix in the range [start, idx1).
        for (j=AR.start; j<idx1; j=AR.succ[j]) {
          f_weighted(&D_(j, idx2), D_(j, idx1) );
          if (n_nghbr[j] == idx1)
            n_nghbr[j] = idx2;
        }
        // Update the distance matrix in the range (idx1, idx2).
        for (; j<idx2; j=AR.succ[j]) {
          f_weighted(&D_(j, idx2), D_(idx1, j) );
          if (D_(j, idx2)<mindist[j]) {
            minidx.update_leq(j, D_(j, idx2));
            n_nghbr[j] = idx2;
          }
        }
        // Update the distance matrix in the range (idx2, N).
        if (idx2<N_1) {
          n_nghbr[idx2] = j = AR.succ[idx2]; // exists, maximally N-1
          f_weighted(&D_(idx2, j), D_(idx1, j) );
          min = D_(idx2,j);
          for (j=AR.succ[j]; j<N; j=AR.succ[j]) {
            f_weighted(&D_(idx2, j), D_(idx1, j) );
            if (D_(idx2,j)<min) {
              min = D_(idx2,j);
              n_nghbr[idx2] = j;
            }
          }
          minidx.update(idx2, min);
        }
        break;
      }

      case METHOD_METR_WARD: {
        /*
        Ward linkage.

        Shorter and longer distances can occur, not smaller than min(d1,d2)
        but maybe bigger than max(d1,d2).
        */
        // Update the distance matrix in the range [start, idx1).
        for (j=AR.start; j<idx1; j=AR.succ[j]) {
          f_ward(&D_(j, idx2), D_(j, idx1), mindist[idx1],
                 size1, size2, static_cast<t_float>(members[j]) );
          if (n_nghbr[j] == idx1)
            n_nghbr[j] = idx2;
        }
        // Update the distance matrix in the range (idx1, idx2).
        for (; j<idx2; j=AR.succ[j]) {
          f_ward(&D_(j, idx2), D_(idx1, j), mindist[idx1], size1, size2,
                 static_cast<t_float>(members[j]) );
          if (D_(j, idx2)<mindist[j]) {
            minidx.update_leq(j, D_(j, idx2));
            n_nghbr[j] = idx2;
          }
        }
        // Update the distance matrix in the range (idx2, N).
        if (idx2<N_1) {
          n_nghbr[idx2] = j = AR.succ[idx2]; // exists, maximally N-1
          f_ward(&D_(idx2, j), D_(idx1, j), mindist[idx1],
                 size1, size2, static_cast<t_float>(members[j]) );
          min = D_(idx2,j);
          for (j=AR.succ[j]; j<N; j=AR.succ[j]) {
            f_ward(&D_(idx2, j), D_(idx1, j), mindist[idx1],
                   size1, size2, static_cast<t_float>(members[j]) );
            if (D_(idx2,j)<min) {
              min = D_(idx2,j);
              n_nghbr[idx2] = j;
            }
          }
          minidx.update(idx2, min);
        }
        break;
      }

      case METHOD_METR_CENTROID: {
        /*
        Centroid linkage.

        Shorter and longer distances can occur, not bigger than max(d1,d2)
        but maybe smaller than min(d1,d2).
        */
        // Update the distance matrix in the range [start, idx1).
        t_float s = size1/(size1+size2);
        t_float t = size2/(size1+size2);
        for (j=AR.start; j<idx1; j=AR.succ[j]) {
          f_centroid(&D_(j, idx2), D_(j, idx1), mindist[idx1], s, t);
          if (D_(j, idx2)<mindist[j]) {
            minidx.update_leq(j, D_(j, idx2));
            n_nghbr[j] = idx2;
          }
          else if (n_nghbr[j] == idx1)
            n_nghbr[j] = idx2;
        }
        // Update the distance matrix in the range (idx1, idx2).
        for (; j<idx2; j=AR.succ[j]) {
          f_centroid(&D_(j, idx2), D_(idx1, j), mindist[idx1], s, t);
          if (D_(j, idx2)<mindist[j]) {
            minidx.update_leq(j, D_(j, idx2));
            n_nghbr[j] = idx2;
          }
        }
        // Update the distance matrix in the range (idx2, N).
        if (idx2<N_1) {
          n_nghbr[idx2] = j = AR.succ[idx2]; // exists, maximally N-1
          f_centroid(&D_(idx2, j), D_(idx1, j), mindist[idx1], s, t);
          min = D_(idx2,j);
          for (j=AR.succ[j]; j<N; j=AR.succ[j]) {
            f_centroid(&D_(idx2, j), D_(idx1, j), mindist[idx1], s, t);
            if (D_(idx2,j)<min) {
              min = D_(idx2,j);
              n_nghbr[idx2] = j;
            }
          }
          minidx.update(idx2, min);
        }
        break;
      }

      case METHOD_METR_MEDIAN: {
        /*
        Median linkage.

        Shorter and longer distances can occur, not bigger than max(d1,d2)
        but maybe smaller than min(d1,d2).
        */
        // Update the distance matrix in the range [start, idx1).
        for (j=AR.start; j<idx1; j=AR.succ[j]) {
          f_median(&D_(j, idx2), D_(j, idx1), mindist[idx1] );
          if (D_(j, idx2)<mindist[j]) {
            minidx.update_leq(j, D_(j, idx2));
            n_nghbr[j] = idx2;
          }
          else if (n_nghbr[j] == idx1)
            n_nghbr[j] = idx2;
        }
        // Update the distance matrix in the range (idx1, idx2).
        for (; j<idx2; j=AR.succ[j]) {
          f_median(&D_(j, idx2), D_(idx1, j), mindist[idx1] );
          if (D_(j, idx2)<mindist[j]) {
            minidx.update_leq(j, D_(j, idx2));
            n_nghbr[j] = idx2;
          }
        }
        // Update the distance matrix in the range (idx2, N).
        if (idx2<N_1) {
          n_nghbr[idx2] = j = AR.succ[idx2]; // exists, maximally N-1
          f_median(&D_(idx2, j), D_(idx1, j), mindist[idx1] );
          min = D_(idx2,j);
          for (j=AR.succ[j]; j<N; j=AR.succ[j]) {
            f_median(&D_(idx2, j), D_(idx1, j), mindist[idx1] );
            if (D_(idx2,j)<min) {
              min = D_(idx2,j);
              n_nghbr[idx2] = j;
            }
          }
          minidx.update(idx2, min);
        }
        break;
      }
    }
  }
}
