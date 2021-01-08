
#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "kdtree/kdtree.h"

template <typename T>
KDTree<T>::KDTree() :
    data_(nullptr), n_points_(0), implicit_idx_tree_(nullptr),
    visited_(0), dimension_(0), copied_(false) {}

// a helper function that calculates whether
// there is an integer n such that num_nodes = 2^n - 1.
// if there is, we can build a complete binary tree.
static bool canBuildCompleteBinaryTree(size_t num_nodes) {
  num_nodes = num_nodes + 1;
  while (num_nodes > 1) {
    if ((num_nodes % 2)!=0) {
      return false;
    }
    num_nodes = num_nodes / 2;
  }
  return true;
}

static size_t calcLeftRightBottomNodeDifference(size_t num_nodes) {
  size_t num_bottom_nodes = num_nodes;
  size_t max_nodes_in_curr_level = 1;
  while (num_bottom_nodes >= max_nodes_in_curr_level) {
    num_bottom_nodes = num_bottom_nodes - max_nodes_in_curr_level;
    max_nodes_in_curr_level = max_nodes_in_curr_level * 2;
  }
  size_t left_num_bottom_nodes = 0;
  size_t right_num_bottom_nodes = 0;
  if (num_bottom_nodes > max_nodes_in_curr_level / 2) {
    left_num_bottom_nodes = max_nodes_in_curr_level / 2;
    right_num_bottom_nodes = num_bottom_nodes - max_nodes_in_curr_level / 2;
  } else {
    left_num_bottom_nodes = num_bottom_nodes;
    right_num_bottom_nodes = 0;
  }
  return left_num_bottom_nodes - right_num_bottom_nodes;
}

// a helper function that qsorts an array until
// array[:idx] <= array[idx] <= array[idx+1:]
// for a given integer idx
template<typename T>
static void qsortForThisIndex(const T* data, size_t* idxarr, int dim, size_t size, int target_dim, size_t target_idx) {
  assert(target_idx>=0);
  assert(target_idx<size);
  if (size==1) {
    return;
  }
  size_t tmp;
  const T pivot = data[idxarr[0]*dim+target_dim];
  size_t less = 1;
  size_t more = size-1;
  while (less <= more) {
    const T curr = data[idxarr[less]*dim+target_dim];
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
    qsortForThisIndex<T>(data, idxarr+pivot_pos+1, dim, size-pivot_pos-1, target_dim, target_idx-pivot_pos-1);
  } else if (pivot_pos > target_idx) {
    qsortForThisIndex<T>(data, idxarr, dim, pivot_pos, target_dim, target_idx);
  }
}

// a helper function for the helper function that builds an implicit KDTree
template <typename T>
static void buildImplicitKDTreeHelper(const T* data, size_t* idxarr1, size_t* idxarr2,
                                      size_t start, size_t end, int dim, int curr_dim, size_t idx) {

  if (start > end){
    return;
  }

  if (start==end) {
    idxarr2[idx] = idxarr1[end];
    return;
  }

  size_t size = end - start + 1;

//   size_t mid = start + (end - start + 1) / 2;
//   while (!canBuildCompleteBinaryTree(mid-start)&&!canBuildCompleteBinaryTree(end-mid)) {
//     mid++;
//   }

  // size_t num_bottom_nodes = (end - start + 1);
  // size_t max_nodes_in_curr_level = 1;
  // while (num_bottom_nodes >= max_nodes_in_curr_level) {
  //   num_bottom_nodes = num_bottom_nodes - max_nodes_in_curr_level;
  //   max_nodes_in_curr_level = max_nodes_in_curr_level * 2;
  // }
  // size_t left_num_bottom_nodes = 0;
  // size_t right_num_bottom_nodes = 0;
  // if (num_bottom_nodes > max_nodes_in_curr_level / 2) {
  //   left_num_bottom_nodes = max_nodes_in_curr_level / 2;
  //   right_num_bottom_nodes = num_bottom_nodes - max_nodes_in_curr_level / 2;
  // } else {
  //   left_num_bottom_nodes = num_bottom_nodes;
  //   right_num_bottom_nodes = 0;
  // }

  size_t diff = calcLeftRightBottomNodeDifference(size);
  size_t mid = start + (size + diff) / 2;

  qsortForThisIndex<T>(data, idxarr1+start, dim, end-start+1, curr_dim, mid-start);

  idxarr2[idx] = idxarr1[mid];
  curr_dim = (curr_dim + 1) % dim;
  buildImplicitKDTreeHelper<T>(data, idxarr1, idxarr2, start, mid-1, dim, curr_dim, idx*2+1);
  buildImplicitKDTreeHelper<T>(data, idxarr1, idxarr2, mid+1, end, dim, curr_dim, idx*2+2);
}

// a helper function that build an implicit KDTree
template <typename T>
static size_t* buildImplicitKDTree(const T* data, int dim, size_t num) {
  size_t* idxarr1 = (size_t*)malloc(sizeof(size_t)*num);
  size_t* idxarr2 = (size_t*)malloc(sizeof(size_t)*num);
  for (size_t i = 0; i < num; i++) {
    idxarr1[i] = i;
  }
  size_t start = 0;
  size_t end = num - 1;
  buildImplicitKDTreeHelper<T>(data,idxarr1,idxarr2,start,end,dim,0,0);
  free(idxarr1);
  return idxarr2;
}

// KDTree constructor with a raw pointer
template <typename T>
KDTree<T>::KDTree(const T* data, int dim, size_t num, bool copy /* = false */) {
    implicit_idx_tree_ = nullptr;
    this->assign(data, dim, num, copy);
}

// KDTree constructor with 1d stdvector
template <typename T>
KDTree<T>::KDTree(const std::vector<T>& data, int dim, bool copy /* false */) {
    implicit_idx_tree_ = nullptr;
    this->assign((const T*)data.data(), dim, data.size() / dim, copy);
}

// destructor
template <typename T>
KDTree<T>::~KDTree() {
  if ((copied_)&&(data_)) { free((T*)data_); }
  if (implicit_idx_tree_) { free(implicit_idx_tree_); }
  // n_points_ = 0;
  // visited_ = 0;
}

// builds a KDTree for a new data
template <typename T>
void KDTree<T>::assign(const T* data, int dim, size_t num, bool copy /* = false */) {
  if ((copied_)&&(data_)) {
      free((T*)data_);
  }
  if (!copy) {
    data_ = data;
    copied_ = false;
  } else {
    data_ = (T*)malloc(sizeof(T)*num);
    memcpy((T*)data_, data, sizeof(T)*num);
    copied_ = true;
  }
  dimension_ = dim;
  n_points_ = num;
  visited_ = 0;
  if (implicit_idx_tree_) { free(implicit_idx_tree_); }
  implicit_idx_tree_ = buildImplicitKDTree<T>(data, dim, num);

  printf("kd tree successfully built.\n");
}

template <typename T>
void KDTree<T>::assign(const std::vector<T>& data, int dim, bool copy  /* = false */) {
  this->assign((const T*)data.data(), dim, data.size() / dim, copy);
}

template <typename T>
bool KDTree<T>::isCopied() {
  return copied_;
}

#include "kdtree.cpp"
