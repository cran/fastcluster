/*
  fastcluster: Fast hierarchical clustering routines for R and Python

  Copyright © 2011 Daniel Müllner
  <http://math.stanford.edu/~muellner>
*/
#define __STDC_LIMIT_MACROS
//#include <stdint.h>

#include <Python.h>
#include <numpy/arrayobject.h>

#include "fastcluster.cpp"

class linkage_output {
private:
  t_float *Z;
  t_index pos;

public:
  linkage_output(t_float * Z) {
    this->Z = Z;
    pos = 0;
  }

  void append(t_index node1, t_index node2, t_float dist, t_float size) {
    if (node1<node2) {
      Z[pos++] = static_cast<t_float>(node1);
      Z[pos++] = static_cast<t_float>(node2);
    }
    else {
      Z[pos++] = static_cast<t_float>(node2);
      Z[pos++] = static_cast<t_float>(node1);
    }
    Z[pos++] = dist;
    Z[pos++] = size;
  }
};

// The size of a node is either 1 (a single point) or is looked up from
// one of the clusters.
#define size_(r_) ( ((r_<N) ? 1 : Z_(r_-N,3)) )

template <bool sorted>
static void generate_SciPy_dendrogram(t_float * const Z, cluster_result & Z2, const t_index N) {
  // The array "nodes" is a union-find data structure for the cluster
  // identites (only needed for unsorted cluster_result input).
  union_find nodes;
  if (!sorted) {
    std::stable_sort(Z2[0], Z2[N-1]);
    nodes.init(N);
  }

  linkage_output output(Z);
  t_index node1, node2;

  for (t_index i=0; i<N-1; i++) {
    // Get two data points whose clusters are merged in step i.
    if (sorted) {
      node1 = Z2[i]->node1;
      node2 = Z2[i]->node2;
    }
    else {
      // Find the cluster identifiers for these points.
      node1 = nodes.Find(Z2[i]->node1);
      node2 = nodes.Find(Z2[i]->node2);
      // Merge the nodes in the union-find data structure by making them
      // children of a new node.
      nodes.Union(node1, node2);
    }
    output.append(node1, node2, Z2[i]->dist, size_(node1)+size_(node2));
  }
}

static void square(t_index N, t_float * const D) {
  for (long_index i=0; i<static_cast<long_index>(N)*(N-1)/2; i++)
    D[i] *= D[i];
}

/*
  Python interface code
*/
static PyObject * linkage_wrap(PyObject *self, PyObject *args);

// List the C++ methods that this extension provides.
static PyMethodDef _fastclusterWrapMethods[] = {
  {"linkage_wrap", linkage_wrap, METH_VARARGS},
  {NULL, NULL, 0, NULL}     /* Sentinel - marks the end of this structure */
};

// Tell Python about these methods.
PyMODINIT_FUNC init_fastcluster(void)  {
  (void) Py_InitModule("_fastcluster", _fastclusterWrapMethods);
  import_array();  // Must be present for NumPy. Called first after above line.
}

// This is the interface to Python.
static PyObject *linkage_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *D, *Z;
  long int N = 0;
  unsigned char method;
  try{
    // Parse the input arguments
    if (!PyArg_ParseTuple(args, "lO!O!b",
                          &N,                // signed long integer
                          &PyArray_Type, &D, // NumPy array
                          &PyArray_Type, &Z, // NumPy array
                          &method)) {        // unsigned char
      return NULL; // Error if the arguments have the wrong type.
    }
    if (N < 1 ) {
      // N must be at least 1.
      PyErr_SetString(PyExc_ValueError,
                      "At least one element is needed for clustering.");
      return NULL;
    }

    // (1)
    // The biggest index used below is 4*(N-2)+3, as an index to Z. This must fit
    // into the data type used for indices.
    // (2)
    // The largest representable integer, without loss of precision, by a floating
    // point number of type t_float is 2^T_FLOAT_MANT_DIG. Here, we make sure that
    // all cluster labels from 0 to 2N-2 in the output can be accurately represented
    // by a floating point number.
    if (N > MAX_INDEX/4 ||
        (N-1)>>(T_FLOAT_MANT_DIG-1) > 0) {
      PyErr_SetString(PyExc_ValueError,
                      "Data is too big, index overflow.");
      return NULL;
    }

    t_float * const D_ = reinterpret_cast<t_float *>(D->data);
    t_float * const Z_ = reinterpret_cast<t_float *>(Z->data);

    cluster_result Z2(N-1);

    auto_array_ptr<t_index> members;
    if (method==METHOD_METR_AVERAGE ||
        method==METHOD_METR_WARD ||
        method==METHOD_METR_CENTROID) {
      members.init(N, 1);
    }

    if (method==METHOD_METR_WARD ||
        method==METHOD_METR_CENTROID ||
        method==METHOD_METR_MEDIAN) {
      square(N, D_);
    }

    switch (method) {
    case METHOD_METR_SINGLE:
      single_linkage(N, D_, Z2);
      break;
    case METHOD_METR_COMPLETE:
      NN_chain_linkage<METHOD_METR_COMPLETE, t_index>(N, D_, NULL, Z2);
      break;
    case METHOD_METR_AVERAGE:
      NN_chain_linkage<METHOD_METR_AVERAGE, t_index>(N, D_, members, Z2);
      break;
    case METHOD_METR_WEIGHTED:
      NN_chain_linkage<METHOD_METR_WEIGHTED, t_index>(N, D_, NULL, Z2);
      break;
    case METHOD_METR_WARD:
      NN_chain_linkage<METHOD_METR_WARD, t_index>(N, D_, members, Z2);
      break;
    case METHOD_METR_CENTROID:
      generic_linkage<METHOD_METR_CENTROID, t_index>(N, D_, members, Z2);
      break;
    case METHOD_METR_MEDIAN:
      generic_linkage<METHOD_METR_MEDIAN, t_index>(N, D_, NULL, Z2);
      break;
    default:
      PyErr_SetString(PyExc_IndexError, "Invalid method index.");
      return NULL;
    }

    if (method==METHOD_METR_AVERAGE ||
        method==METHOD_METR_WARD ||
        method==METHOD_METR_CENTROID) {
      members.free();
    }

    if (method==METHOD_METR_WARD ||
        method==METHOD_METR_CENTROID ||
        method==METHOD_METR_MEDIAN) {
      Z2.sqrt();
    }

    if (method==METHOD_METR_CENTROID ||
        method==METHOD_METR_MEDIAN) {
      generate_SciPy_dendrogram<true>(Z_, Z2, N);
    }
    else {
      generate_SciPy_dendrogram<false>(Z_, Z2, N);
    }

    Py_INCREF(Py_None); // Return None on success.
  } // try
  catch (std::bad_alloc&) {
    return PyErr_NoMemory();
  }
  catch(std::exception& e){
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  catch(...){
    PyErr_SetString(PyExc_StandardError,
                    "C++ exception (unknown reason). Please send a bug report.");
    return NULL;
  }

  return Py_None; // see Py_INCREF(Py_None) above
}
