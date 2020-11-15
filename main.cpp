#include <iostream>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include "geometry3d.hpp"
#include "utils.hpp"
#include <omp.h>

using namespace std;

int main(int argc, char** argv) {
  if (argc < 2) {
    printf("Not Enough Arguments\n");
    exit(1);
  }

  char* path = argv[1];
  Geometry3D pcd(path);
  printf("file name: %s\n", pcd.path);
  printf("number of points: %ld\n", pcd.n_points);
  float* points = pcd.points;

  pcd.build_KDTree();
  int k = 20;
  // float arr[3] = {-72.202599,-51.793232, 525.919250};
  float arr[3] = {-21.996,-25.391,6.480};
  size_t* k_nearest = new size_t[k];
  float* k_distances = new float[k];

  int visited = pcd.find_k_nearest(arr,k,k_nearest,k_distances);
  printf("visited: %d\n",visited);
  printf("%d nearest points\n",k);
  for (int i = 0; i < k; i++) {
    size_t idx = k_nearest[i];
    float* x = pcd.points + idx * 3;
    printf("%8ld: <%.3f,%.3f,%.3f>, ",idx,x[0],x[1],x[2]);
    printf("d = %.5f\n",k_distances[i]);
  }

  delete[] k_nearest;
  delete[] k_distances;
  return 0;
}
