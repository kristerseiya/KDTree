
#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "kdtree/kdtree.h"

template <typename T>
KDTree<T>::KDTree() :
    data_(nullptr), n_points_(0), implicit_idx_tree_(nullptr),
    visited_(0), dimension_(0), copied_(false), leaf_size_(1),
    leaf_starts_at_(0) {}

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

static int calcLeftRightBottomNodeDifferenceLeafOne(size_t num_nodes) {
  int num_bottom_nodes = num_nodes;
  int max_nodes_in_curr_level = 1;
  while (num_bottom_nodes >= max_nodes_in_curr_level) {
    num_bottom_nodes = num_bottom_nodes - max_nodes_in_curr_level;
    max_nodes_in_curr_level = max_nodes_in_curr_level * 2;
  }
  int left_num_bottom_nodes = 0;
  int right_num_bottom_nodes = 0;
  if (num_bottom_nodes > max_nodes_in_curr_level / 2) {
    left_num_bottom_nodes = max_nodes_in_curr_level / 2;
    right_num_bottom_nodes = num_bottom_nodes - max_nodes_in_curr_level / 2;
  } else {
    left_num_bottom_nodes = num_bottom_nodes;
    right_num_bottom_nodes = 0;
  }
  return left_num_bottom_nodes - right_num_bottom_nodes;
}

static size_t computeLeafNodeStartIndex(size_t num_nodes, size_t leaf_size) {
    // if (leaf_size == 1) {
    //   num_nodes = num_nodes + 1;
    // }
    size_t num_leaf_nodes = (num_nodes + 1) / (1 + leaf_size);
    return num_leaf_nodes - 1;
}

static int calcLeftRightBottomNodeDifference(size_t num_nodes, size_t leaf_size) {

  // if (leaf_size == 1) {
  //   return calcLeftRightBottomNodeDifferenceLeafOne(num_nodes);
  // }

  int num_leaf_nodes = (num_nodes + 1) / (1 + leaf_size);
  int num_bottom_nodes = 2 * num_leaf_nodes - 1;
  int max_nodes_in_curr_level = 1;

  while (num_bottom_nodes >= max_nodes_in_curr_level) {
    num_bottom_nodes = num_bottom_nodes - max_nodes_in_curr_level;
    max_nodes_in_curr_level = max_nodes_in_curr_level * 2;
  }

  int left_num_bottom_nodes = 0;
  int right_num_bottom_nodes = 0;
  if (num_bottom_nodes > max_nodes_in_curr_level / 2) {
    left_num_bottom_nodes = max_nodes_in_curr_level / 2;
    right_num_bottom_nodes = num_bottom_nodes - max_nodes_in_curr_level / 2;
  } else {
    left_num_bottom_nodes = num_bottom_nodes;
    right_num_bottom_nodes = 0;
  }

  int left_extra = left_num_bottom_nodes * leaf_size
              + (max_nodes_in_curr_level / 2 - left_num_bottom_nodes) / 2 * (leaf_size - 1);
  int right_extra = right_num_bottom_nodes * leaf_size
              + (max_nodes_in_curr_level / 2 - right_num_bottom_nodes) / 2 * (leaf_size - 1);
  int extra = num_nodes - (num_leaf_nodes - 1) - num_leaf_nodes * leaf_size;

  return (right_num_bottom_nodes == 0) ? left_extra - right_extra + extra : left_extra - right_extra - extra;
}

static int calcSplitLeafOne(size_t num_nodes) {
  int num_bottom_nodes = num_nodes;
  int max_nodes_in_curr_level = 1;
  while (num_bottom_nodes >= max_nodes_in_curr_level) {
    num_bottom_nodes = num_bottom_nodes - max_nodes_in_curr_level;
    max_nodes_in_curr_level = max_nodes_in_curr_level * 2;
  }
  int left_num_bottom_nodes = 0;
  int right_num_bottom_nodes = 0;
  if (num_bottom_nodes > max_nodes_in_curr_level / 2) {
    left_num_bottom_nodes = max_nodes_in_curr_level / 2;
    right_num_bottom_nodes = num_bottom_nodes - max_nodes_in_curr_level / 2;
  } else {
    left_num_bottom_nodes = num_bottom_nodes;
    right_num_bottom_nodes = 0;
  }
  return left_num_bottom_nodes + (max_nodes_in_curr_level - 2) / 2;
}

static size_t calcSplit(size_t num_nodes, size_t leaf_size) {

  // if (leaf_size == 1) {
  //     return calcSplitLeafOne(num_nodes);
  // }

  size_t num_leaf_nodes = (num_nodes + 1) / (1 + leaf_size);
  size_t num_bottom_nodes = 2 * num_leaf_nodes - 1;
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

  size_t left_extra = left_num_bottom_nodes * leaf_size
              + (max_nodes_in_curr_level / 2 - left_num_bottom_nodes) / 2 * (leaf_size - 1);
  size_t right_extra = right_num_bottom_nodes * leaf_size
              + (max_nodes_in_curr_level / 2 - right_num_bottom_nodes) / 2 * (leaf_size - 1);
  size_t extra = num_nodes - (num_leaf_nodes - 1) - num_leaf_nodes * leaf_size;

  return ((right_num_bottom_nodes == 0)&&(left_num_bottom_nodes!=0)) ?
  left_extra + extra + (max_nodes_in_curr_level - 2) / 2 : left_extra + (max_nodes_in_curr_level - 2) / 2;
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
                                      size_t start, size_t end, int dim, int curr_dim,
                                      size_t idx, const size_t leaf_size, const size_t leaf_node_starts) {

  if (start > end) {
    return;
  }

  size_t size = end - start + 1;

  if (idx >= leaf_node_starts) {
    assert(size >= leaf_size);
    // if (size >= 2 * leaf_size) {
    //   printf("%lu\n", size);
    //   printf("%lu\n", idx);
    //   printf("%lu\n", leaf_node_starts);
    // }
    assert(size <= 2 * leaf_size);
    idx = leaf_node_starts + (idx - leaf_node_starts) * leaf_size;
    for (int i = 0; i < size; i++) {
      idxarr2[idx+i] = idxarr1[start+i];
    }
    return;
  }
  // if (start == end) {
  //   idxarr2[idx] = idxarr1[end];
  //   return;
  // }

//   size_t mid = start + (end - start + 1) / 2;
//   while (!canBuildCompleteBinaryTree(mid-start)&&!canBuildCompleteBinaryTree(end-mid)) {
//     mid++;
//   }

  size_t mid = start + calcSplit(size, leaf_size);
  // int diff = calcLeftRightBottomNodeDifference(size, leaf_size);
  // size_t mid = start + (size + diff) / 2;
  // size_t mid = start;
  // while ((mid - start) - (end - mid) != diff) {
  //   mid++;
  // }

  qsortForThisIndex<T>(data, idxarr1+start, dim, end-start+1, curr_dim, mid-start);

  idxarr2[idx] = idxarr1[mid];
  curr_dim = (curr_dim + 1) % dim;
  buildImplicitKDTreeHelper<T>(data, idxarr1, idxarr2, start, mid-1, dim, curr_dim, idx*2+1, leaf_size, leaf_node_starts);
  buildImplicitKDTreeHelper<T>(data, idxarr1, idxarr2, mid+1, end, dim, curr_dim, idx*2+2, leaf_size, leaf_node_starts);
}

// a helper function that build an implicit KDTree
template <typename T>
static size_t* buildImplicitKDTree(const T* data, int dim, size_t num, const size_t leaf_size, const size_t leaf_node_starts) {
  size_t* idxarr1 = (size_t*)malloc(sizeof(size_t)*num);
  size_t* idxarr2 = (size_t*)malloc(sizeof(size_t)*num);
  for (size_t i = 0; i < num; i++) {
    idxarr1[i] = i;
  }
  size_t start = 0;
  size_t end = num - 1;
  buildImplicitKDTreeHelper<T>(data, idxarr1, idxarr2, start, end, dim, 0, 0, leaf_size, leaf_node_starts);
  free(idxarr1);
  return idxarr2;
}

// KDTree constructor with a raw pointer
template <typename T>
KDTree<T>::KDTree(const T* data, int dim, size_t num, int leaf_size /* = 1 */, bool copy /* = false */) {
    implicit_idx_tree_ = nullptr;
    this->assign(data, dim, num, leaf_size, copy);
}

// KDTree constructor with 1d stdvector
template <typename T>
KDTree<T>::KDTree(const std::vector<T>& data, int dim, int leaf_size /* = 1 */, bool copy /* false */) {
    implicit_idx_tree_ = nullptr;
    this->assign((const T*)data.data(), dim, data.size() / dim, leaf_size, copy);
}

// destructor
template <typename T>
KDTree<T>::~KDTree() {
  if ((copied_)&&(data_)) { free((T*)data_); }
  if (implicit_idx_tree_) { free(implicit_idx_tree_); }
}

// builds a KDTree for a new data
template <typename T>
void KDTree<T>::assign(const T* data, int dim, size_t num, int leaf_size /* = 1 */, bool copy /* = false */) {
  if (num < 1) {
    throw std::runtime_error("data size must be bigger than 0");
  }
  if (dim < 1) {
    throw std::runtime_error("dimenstion must be bigger than 0");
  }
  if (leaf_size < 1) {
    throw std::runtime_error("leaf size must be integer bigger than 0");
  }

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
  leaf_size_ = leaf_size;
  leaf_starts_at_ = computeLeafNodeStartIndex(n_points_, leaf_size_);

  if (implicit_idx_tree_) { free(implicit_idx_tree_); }
  implicit_idx_tree_ = buildImplicitKDTree<T>(data, dim, num, leaf_size_, leaf_starts_at_);

  printf("kd tree successfully built.\n");
}

template <typename T>
void KDTree<T>::assign(const std::vector<T>& data, int dim, int leaf_size /* = 1 */, bool copy  /* = false */) {
  this->assign((const T*)data.data(), dim, data.size() / dim, leaf_size, copy);
}

template <typename T>
bool KDTree<T>::isCopied() {
  return copied_;
}

#include "kdtree.cpp"
