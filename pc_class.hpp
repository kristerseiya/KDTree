#ifndef PC_CLASS_H
#define PC_CLASS_H

#include <stdlib.h>
#include <stdbool.h>
#include "kdtree.hpp"

using namespace std;
typedef unsigned char uint8;

class PointCloud {
public:
  char* path;
  size_t num_points;
  size_t width;
  size_t height;
  bool is_kdtree_built;
  bool has_normals;
  PointCloud();
  PointCloud(char* path/*, bool preprocess_nan = true */);
  ~PointCloud();
  void build_KDTree();
  int find_k_nearest(float* query, int k, size_t* k_nearest, float* k_distances);

public:
  float* points;
  float* normals;
  uint8* colors;
  // uint8* mask;

private:
  // bool preprocess_nan;
  KDTree kdtree;
  void read_ply(char* path);
  void read_xyzm(char* paht);
};

#endif
