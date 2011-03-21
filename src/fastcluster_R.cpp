/*
  fastcluster: Fast hierarchical clustering routines for R and Python

  Copyright © 2011 Daniel Müllner
  <http://math.stanford.edu/~muellner>
*/
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "fastcluster.cpp"

struct node_pos {
  int node;
  t_index pos;
};

void order_nodes(const int N, const int * const merge, const t_index * const node_size, int * const order) {
  auto_array_ptr<node_pos> queue(N/2);

  int parent;
  int child;
  t_index pos = 0;

  queue[0].node = N-2;
  queue[0].pos = 0;
  t_index idx = 1;

  do {
    idx --;
    parent = queue[idx].node;
    pos = queue[idx].pos;

    // First child
    child = merge[parent];
    if (child<0) {
      order[pos] = -child;
      pos++;
    }
    else {
      queue[idx].node = child-1;
      queue[idx].pos = pos;
      idx++;
      pos += node_size[child-1];
    }
    // Second child
    child = merge[parent+N-1];
    if (child<0) {
      order[pos] = -child;
    }
    else {
      queue[idx].node = child-1;
      queue[idx].pos = pos;
      idx++;
    }
  } while (idx>0);
}

#define size_(r_) ( ((r_<N) ? 1 : node_size[r_-N]) )

template <bool sorted>
void generate_R_dendrogram(int * merge, double * const height, int * const order, cluster_result & Z2, const int N) {
  // The array "nodes" is a union-find data structure for the cluster
  // identites (only needed for unsorted cluster_result input).
  union_find nodes;

  if (!sorted) {
    std::stable_sort(Z2[0], Z2[N-1]);
    nodes.init(N);
  }

  t_index node1, node2;
  auto_array_ptr<t_index> node_size(N-1);

  for (t_index i=0; i<N-1; i++) {
    // Get two data points whose clusters are merged in step i.
    // Find the cluster identifiers for these points.
    if (sorted) {
      node1 = Z2[i]->node1;
      node2 = Z2[i]->node2;
    }
    else {
      node1 = nodes.Find(Z2[i]->node1);
      node2 = nodes.Find(Z2[i]->node2);
      // Merge the nodes in the union-find data structure by making them
      // children of a new node.
      nodes.Union(node1, node2);
    }
    // Sort the nodes in the output array.
    if (node1>node2) {
      t_index tmp = node1;
      node1 = node2;
      node2 = tmp;
    }
    merge[i]     = (node1<N) ? -static_cast<int>(node1)-1
                              : static_cast<int>(node1)-N+1;
    merge[i+N-1] = (node2<N) ? -static_cast<int>(node2)-1
                              : static_cast<int>(node2)-N+1;
    height[i] = Z2[i]->dist;
    node_size[i] = size_(node1) + size_(node2);
  }

  order_nodes(N, merge, node_size, order);
}

/*
  R interface code
*/

extern "C" {
  SEXP fastcluster(SEXP N_, SEXP method_, SEXP D_, SEXP members_) {
    SEXP r = NULL; // return value

    try{
      /*
        Input checks
      */
      // Parameter N: number of data points
      PROTECT(N_);
      if (!IS_INTEGER(N_) || LENGTH(N_)!=1)
        Rf_error("'N' must be a single integer.");
      const int N = *INTEGER_POINTER(N_);
      if (N<2)
        Rf_error("N must be at least 2.");
      const long_index NN = static_cast<long_index>(N)*(N-1)/2;
      UNPROTECT(1); // N_

      // Parameter method: dissimilarity index update method
      PROTECT(method_);
      if (!IS_INTEGER(method_) || LENGTH(method_)!=1) {
        Rf_error("'method' must be a single integer.");
      }
      const int method = *INTEGER_POINTER(method_) - 1; // index-0 based;
      if (method<METHOD_METR_SINGLE || method>METHOD_METR_MEDIAN) {
        Rf_error("Invalid method index.");
      }
      UNPROTECT(1); // method_

      // Parameter members: number of members in each node
      auto_array_ptr<int> members;
      if (method==METHOD_METR_AVERAGE ||
          method==METHOD_METR_WARD ||
          method==METHOD_METR_CENTROID) {
        members.init(N);
        if (Rf_isNull(members_)) {
          for (t_index i=0; i<N; i++) members[i] = 1;
        }
        else {
          PROTECT(members_ = AS_INTEGER(members_));
          if (LENGTH(members_)!=N)
            Rf_error("'members' must have length N.");
          const int * const m = INTEGER_POINTER(members_);
          for (t_index i=0; i<N; i++) members[i] = m[i];
          UNPROTECT(1); // members
        }
      }

      // Parameter D_: dissimilarity matrix
      PROTECT(D_ = AS_NUMERIC(D_));
      if (LENGTH(D_)!=NN)
        Rf_error("'D' must have length (N \\choose 2).");
      const double * const D = NUMERIC_POINTER(D_);
      UNPROTECT(1); // D_

      // Make a working copy of the dissimilarity array
      // for all methods except "single".
      auto_array_ptr<double> D__;
      if (method!=METHOD_METR_SINGLE) {
        D__.init(NN);
        for (long_index i=0; i<NN; i++)
          D__[i] = D[i];
      }

      /*
        Clustering step
      */
      cluster_result Z2(N-1);
      switch (method) {
      case METHOD_METR_SINGLE:
        single_linkage(N, D, Z2);
        break;
      case METHOD_METR_COMPLETE:
        NN_chain_linkage<METHOD_METR_COMPLETE, int>(N, D__, NULL, Z2);
        break;
      case METHOD_METR_AVERAGE:
        NN_chain_linkage<METHOD_METR_AVERAGE, int>(N, D__, members, Z2);
        break;
      case METHOD_METR_WEIGHTED:
        NN_chain_linkage<METHOD_METR_WEIGHTED, int>(N, D__, NULL, Z2);
        break;
      case METHOD_METR_WARD:
        NN_chain_linkage<METHOD_METR_WARD, int>(N, D__, members, Z2);
        break;
      case METHOD_METR_CENTROID:
        generic_linkage<METHOD_METR_CENTROID, int>(N, D__, members, Z2);
        break;
      case METHOD_METR_MEDIAN:
        generic_linkage<METHOD_METR_MEDIAN, int>(N, D__, NULL, Z2);
        break;
      }

      D__.free();     // Free the memory now
      members.free(); // (not strictly necessary).

      SEXP m; // return field "merge"
      PROTECT(m = NEW_INTEGER(2*(N-1)));
      int * const merge = INTEGER_POINTER(m);

      SEXP dim; // Specify that m is an (N-1)×2 matrix
      PROTECT(dim = NEW_INTEGER(2));
      INTEGER(dim)[0] = N-1;
      INTEGER(dim)[1] = 2;
      SET_DIM(m, dim);

      SEXP h; // return field "height"
      PROTECT(h = NEW_NUMERIC(N-1));
      double * const height = NUMERIC_POINTER(h);

      SEXP o; // return fiels "order'
      PROTECT(o = NEW_INTEGER(N));
      int * const order = INTEGER_POINTER(o);

      if (method==METHOD_METR_CENTROID ||
          method==METHOD_METR_MEDIAN)
        generate_R_dendrogram<true>(merge, height, order, Z2, N);
      else
        generate_R_dendrogram<false>(merge, height, order, Z2, N);

      SEXP n; // names
      PROTECT(n = NEW_CHARACTER(3));
      SET_STRING_ELT(n, 0, COPY_TO_USER_STRING("merge"));
      SET_STRING_ELT(n, 1, COPY_TO_USER_STRING("height"));
      SET_STRING_ELT(n, 2, COPY_TO_USER_STRING("order"));

      PROTECT(r = NEW_LIST(3)); // field names in the output list
      SET_ELEMENT(r, 0, m);
      SET_ELEMENT(r, 1, h);
      SET_ELEMENT(r, 2, o);
      SET_NAMES(r, n);

      UNPROTECT(6); // m, dim, h, o, r, n
    } // try
    catch (std::bad_alloc&) {
      Rf_error( "Memory overflow.");
    }
    catch(std::exception& e){
      Rf_error( e.what() );
    }
    catch(...){
      Rf_error( "C++ exception (unknown reason)." );
    }

    return r;
  }

  void R_init_fastcluster(DllInfo *info)
  {
    R_CallMethodDef callMethods[]  = {
      {"fastcluster", (DL_FUNC) &fastcluster, 4},
      {NULL, NULL, 0}
    };
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  }

} // extern "C"
