#ifndef KDTREE_H
#define KDTREE_H

#include <stdlib.h>

using namespace std;
typedef unsigned char uint8;

class KDTree {
public:
  KDTree();
  KDTree(float* data, int dim, size_t num);
  ~KDTree();
  size_t num_nodes;
  int dimension;
  size_t visited;
  int find_k_nearest(float* query, int k, size_t* k_nearest_idx, float* k_distances);
  void set(float* data, int dim, size_t num);

private:
  float* data;
  size_t* implicit_tree;
  void find_k_nearest_with_implicit_tree(size_t node, float* query, int k, int curr_dim, size_t* k_nearest, float* k_distances, int* found);
};

#endif
