#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "kdtree.hpp"

using namespace std;

KDTree::KDTree() {
    this->data = NULL;
    this->num_nodes = 0;
    this->implicit_tree = NULL;
}

// static void Merge(float* data, size_t* idxarr, size_t start, size_t mid, size_t end, int dim, int curr_dim) {
//   size_t size = end - start + 1;
//   size_t* tmp = new size_t[size];
//   size_t idx1 = start;
//   size_t idx2 = mid;
//   for (int i = 0; i < size; i++) {
//     assert(idx1<mid||idx2<=end);
//     if (idx1 >= mid) {
//       tmp[i] = idxarr[idx2];
//       idx2++;
//     } else if (idx2 > end) {
//       tmp[i] = idxarr[idx1];
//       idx1++;
//     } else {
//       float curr1 = data[idxarr[idx1]*dim+curr_dim];
//       float curr2 = data[idxarr[idx2]*dim+curr_dim];
//       curr1 = isnan(curr1) ? FLT_MAX : curr1;
//       curr2 = isnan(curr2) ? FLT_MAX : curr2;
//       if ((curr1 <= curr2) && (idx1 < mid)) {
//         tmp[i] = idxarr[idx1];
//         idx1++;
//       } else if (!(curr1 <= curr2) && (idx2 <= end)) {
//         tmp[i] = idxarr[idx2];
//         idx2++;
//       }
//     }
//   }
//   for (int i = 0; i < size; i++) {
//     idxarr[start+i] = tmp[i];
//   }
//   delete[] tmp;
//   return;
// }
//
// static void MergeSort(float* data, size_t* idxarr, size_t start, size_t end, int dim, int curr_dim) {
//   if (start >= end) {
//     return;
//   }
//   size_t mid = start + (end - start + 1) / 2;
//   MergeSort(data, idxarr, start, mid-1, dim, curr_dim);
//   MergeSort(data, idxarr, mid, end, dim, curr_dim);
//   Merge(data, idxarr, start, mid, end, dim, curr_dim);
//   return;
// }

static bool canBuildStrictBinaryTree(size_t num_nodes) {
  num_nodes = num_nodes + 1;
  while (num_nodes > 1) {
    if ((num_nodes % 2)!=0) {
      return false;
    }
    num_nodes = num_nodes / 2;
  }
  return true;
}

static void find_median_at_this_idx(float* data, size_t* idxarr, int dim, size_t size, int target_dim, size_t target_idx) {
  assert(target_idx>=0);
  assert(target_idx<size);
  if (size==1) {
    return;
  }
  size_t tmp;
  float pivot = data[idxarr[0]*dim+target_dim];
  pivot = isnan(pivot) ? FLT_MAX : pivot;
  size_t less = 1;
  size_t more = size-1;
  while (less <= more) {
    float curr = data[idxarr[less]*dim+target_dim];
    curr = isnan(curr) ? FLT_MAX : curr;
    if (curr <= pivot) {
      less++;
    } else {
      tmp = idxarr[more];
      idxarr[more] = idxarr[less];
      idxarr[less] = tmp;
      more--;
    }
  }
  size_t pivot_pos = less - 1;
  assert(pivot_pos<size);
  tmp = idxarr[pivot_pos];
  idxarr[pivot_pos] = idxarr[0];
  idxarr[0] = tmp;
  if (pivot_pos < target_idx) {
    find_median_at_this_idx(data, idxarr+pivot_pos+1, dim, size-pivot_pos-1, target_dim, target_idx-pivot_pos-1);
  } else if (pivot_pos > target_idx) {
    find_median_at_this_idx(data, idxarr, dim, pivot_pos, target_dim, target_idx);
  }
  return;
}

static void build_implicit_KDTree_helper(float* data, size_t* idxarr1, size_t* idxarr2, size_t start, size_t end, int dim, int curr_dim, size_t idx) {

  if (start > end){
    return;
  }
  if (start==end) {
    idxarr2[idx] = idxarr1[end];
    return;
  }

  size_t mid = start + (end - start + 1) / 2;
  while (!canBuildStrictBinaryTree(mid-start)&&!canBuildStrictBinaryTree(end-mid)) {
    mid++;
  }

  // MergeSort(data, idxarr1, start, end, dim, curr_dim);
  find_median_at_this_idx(data, idxarr1+start, dim, end-start+1, curr_dim, mid-start);

  idxarr2[idx] = idxarr1[mid];
  curr_dim = (curr_dim + 1) % dim;
  build_implicit_KDTree_helper(data, idxarr1, idxarr2, start, mid-1, dim, curr_dim, idx*2+1);
  build_implicit_KDTree_helper(data, idxarr1, idxarr2, mid+1, end, dim, curr_dim, idx*2+2);
  return;
}

static size_t* build_implicit_KDTree(float* data, int dim, size_t num) {
  size_t* idxarr1 = new size_t[num];
  size_t* idxarr2 = new size_t[num];
  for (size_t i = 0; i < num; i++) {
    idxarr1[i] = i;
  }
  size_t start = 0;
  size_t end = num - 1;
  build_implicit_KDTree_helper(data,idxarr1,idxarr2,start,end,dim,0,0);
  delete[] idxarr1;
  return idxarr2;
}

KDTree::KDTree(float* data, int dim, size_t num) {
  this->data = data;
  this->visited = 0;
  this->num_nodes = num;
  this->dimension = dim;
  this->implicit_tree = build_implicit_KDTree(data,dim,num);
  printf("kd tree successfully built.\n");
}

KDTree::~KDTree() {
  if (this->implicit_tree) { delete[] this->implicit_tree; }
  this->num_nodes = 0;
  this->visited = 0;
}

void KDTree::set(float* data, int dim, size_t num) {
  this->data = data;
  this->dimension = dim;
  this->num_nodes = num;
  this->visited = 0;
  if (this->implicit_tree) { delete[] this->implicit_tree; }
  this->implicit_tree = build_implicit_KDTree(data,dim,num);
  printf("kd tree successfully built.\n");
}

void KDTree::find_k_nearest_with_implicit_tree(size_t node, float* query, int k, int curr_dim, size_t* k_nearest, float* k_distances, int* found) {
  if (node >= this->num_nodes) {
    return;
  }
  int dim = this->dimension;
  this->visited++;
  float* curr_x = this->data + this->implicit_tree[node] * dim;
  float distance = 0;

  for (int i = 0; i < dim; i++) {
    float a = curr_x[i] - query[i];
    distance += a * a;
  }

  distance = isnan(distance) ? FLT_MAX : distance;

  if (*found < k) {
    k_distances[*found] = distance;
    k_nearest[*found] = this->implicit_tree[node];
    (*found)++;
  } else if (distance < k_distances[k-1]) {
    k_distances[k-1] = distance;
    k_nearest[k-1] = this->implicit_tree[node];
  }

  int idx = (*found) - 1;
  while ((idx>0)&&(distance < k_distances[idx-1])) {
    k_distances[idx] = k_distances[idx-1];
    k_nearest[idx] = k_nearest[idx-1];
    idx--;
  }
  if (idx != (*found)-1) {
    k_distances[idx] = distance;
    k_nearest[idx] = this->implicit_tree[node];
  }

  int next_dim = (curr_dim + 1) % dim;
  float dist2bb;
  if (query[curr_dim] <= curr_x[curr_dim]) {
    KDTree::find_k_nearest_with_implicit_tree(node*2+1, query, k, next_dim, k_nearest, k_distances, found);
    dist2bb = query[curr_dim] - curr_x[curr_dim];
    dist2bb = dist2bb * dist2bb;
    dist2bb = isnan(dist2bb) ? FLT_MAX : dist2bb;
    if (((*found)<k)||(dist2bb < k_distances[*found-1])) {
      KDTree::find_k_nearest_with_implicit_tree(node*2+2, query, k, next_dim, k_nearest, k_distances, found);
    }
  } else {
    KDTree::find_k_nearest_with_implicit_tree(node*2+2, query, k, next_dim, k_nearest, k_distances, found);
    dist2bb = query[curr_dim] - curr_x[curr_dim];
    dist2bb = dist2bb * dist2bb;
    dist2bb = isnan(dist2bb) ? FLT_MAX : dist2bb;
    if (((*found)<k)||(dist2bb < k_distances[*found-1])) {
      KDTree::find_k_nearest_with_implicit_tree(node*2+1, query, k, next_dim, k_nearest, k_distances, found);
    }
  }
  return;
}

int KDTree::find_k_nearest(float* query, int k, size_t* k_nearest_idx, float* k_distances) {
  this->visited = 0;
  for (int i = 0; i < k; i++) {
    k_distances[i] = -1;
  }
  int found = 0;
  KDTree::find_k_nearest_with_implicit_tree(0, query, k, 0, k_nearest_idx, k_distances, &found);
  return this->visited;
}
