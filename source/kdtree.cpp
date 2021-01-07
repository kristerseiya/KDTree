
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

// add new data
template <typename T>
void KDTree<T>::add(const T* data, size_t num) {
  // if (!copied_) {
  //   T* new_data_ = (T*)malloc(sizeof(T)*(n_points_ + num * dimension_));
  //   memcpy(new_data_, data_, sizeof(T)*n_points_);
  //   data_ = new_data_;
  //   copied_ = true;
  // } else {
  //   data_ = (T*)realloc(data_, num_nodes + num * dimension_);
  // }
  // for (size_t i = 0; i < num * dimension_; i++) {
  //   data_[n_points_+i] = data[i];
  // }
  //
  // int curr_leftover = n_points_;
  // int n_points_in_curr_level = 1;
  // while (curr_leftover >= n_points_in_curr_level) {
  //   curr_leftover = curr_leftover - n_points_in_curr_level;
  //   n_points_in_curr_level = n_points_in_curr_level * 2;
  // }
  // int to_right = 0;
  // int to_left = 0;
  // if (curr_leftover > n_points_in_curr_level / 2) {
  //   to_left -= n_points_in_curr_level / 2;
  //   to_right -= curr_leftover - n_points_in_curr_level / 2;
  // } else {
  //   to_left -= curr_leftover;
  //   to_right -= 0;
  // }
  // int new_leftover = curr_leftover + num;
  // while (new_leftover > n_points_in_curr_level) {
  //   new_leftover = new_leftover - n_points_in_curr_level;
  //   n_points_in_curr_level = n_points_in_curr_level * 2;
  // }
  // if (new_leftover > n_points_in_curr_level / 2) {
  //   to_left += n_points_in_curr_level / 2;
  //   to_right += new_leftover - n_points_in_curr_level / 2;
  // } else {
  //   to_left += new_leftover;
  //   to_right += 0;
  // }
  //
  // size_t mid = start + (end - start + 1 + left_num_bottom_nodes - right_num_bottom_nodes) / 2;
  // int median_idx = (1 + num + to_left - to_right) / 2;

}

template <typename T>
bool KDTree<T>::isCopied() {
  return copied_;
}

// a helper function that performs insertion sort for the last element
template<typename T>
static void sortLastElementBasedOnDistance(std::vector<size_t>& idx_arr,
                                           std::vector<T>& distances) {
    int idx = distances.size() - 1;
    T last_distance = distances.back();
    size_t last_idx = idx_arr.back();
    while ((idx>0)&&(last_distance < distances[idx-1])) {
        distances[idx] = distances[idx-1];
        idx_arr[idx] = idx_arr[idx-1];
        idx--;
    }
    if (idx != distances.size() - 1) {
        distances[idx] = last_distance;
        idx_arr[idx] = last_idx;
    }
}

template <typename T>
static void swapLastTwo(std::vector<size_t>& idx_arr,
                        std::vector<T>& distances) {
    if (distances.size() < 2) {
      return;
    }

    if ( distances[distances.size()-2] > distances.back() ) {
        T tmp = distances[distances.size()-2];
        distances[distances.size()-2] = distances.back();
        distances.back() = tmp;
        size_t tmp2 = idx_arr[idx_arr.size()-2];
        idx_arr[idx_arr.size()-2] = idx_arr.back();
        idx_arr.back() = tmp2;
    }
}

template <typename T>
static void swapWithMax(std::vector<size_t>& idx_arr,
                        std::vector<T>& distances) {
    if (distances.size() < 2) {
      return;
    }

    T max_dist = distances[0];
    size_t max_idx = 0;
    for (size_t i = 1; i < distances.size() - 1; i++) {
        if (max_dist < distances[i]) {
          max_dist = distances[i];
          max_idx = i;
        }
    }

    distances[max_idx] = distances.back();
    distances.back() = max_dist;
    size_t tmp = idx_arr[max_idx];
    idx_arr[max_idx] = idx_arr.back();
    idx_arr.back() = tmp;
}

template <typename T>
static void updateKNN(const std::vector<T>& query,
                      const T* x, size_t x_idx,
                      std::vector<size_t>& neighbor_idx,
                      std::vector<T>& distances, int k) {
    // compute L2-distance^2 from the current vector to the query
    T distance = 0;
    for (int i = 0; i < query.size(); i++) {
        T a = x[i] - query[i];
        distance += a * a;
    }

    // update k-neighbors
    if (neighbor_idx.size() < k) {
        distances.push_back(distance);
        neighbor_idx.push_back(x_idx);
        sortLastElementBasedOnDistance<T>(neighbor_idx, distances);
        //swapLastTwo<T>(neighbor_idx, distances);
    } else if (distance < distances.back()) {
        distances.back() = distance;
        neighbor_idx.back() = x_idx;
        sortLastElementBasedOnDistance<T>(neighbor_idx, distances);
        //swapWithMax<T>(neighbor_idx, distances);
    }
}

// a helper function for knn search
template<typename T>
void KDTree<T>::searchKNNWithImplicitTreeOtherBranch(size_t node_idx,
                                          const std::vector<T>& query,
                                          const int k, int curr_dim,
                                          std::vector<size_t>& neighbor_idx,
                                          std::vector<T>& distances,
                                          T* dist2bbarr) {

    // index out of range of the implicit tree
    if (node_idx >= n_points_) {
        return;
    }

    // copy
    T* dist2bbarr_cpy = (T*)malloc(sizeof(T)*dimension_);
    for (int i = 0; i < dimension_; i++) {
        dist2bbarr_cpy[i] = dist2bbarr[i];
    }

    // update number of nodes visited
    visited_++;

    // get current vector of interest
    const T* curr_x = data_ + implicit_idx_tree_[node_idx] * dimension_;

    int next_dim = (curr_dim + 1) % dimension_;

    if (query[curr_dim] <= curr_x[curr_dim]) {
        // advance to left in kdtree
        searchKNNWithImplicitTreeOtherBranch(node_idx * 2 + 1, query, k, next_dim, neighbor_idx, distances, dist2bbarr_cpy);
        // compute L2-distance^2 from the boundary to the query
        // T dist2bb = query[curr_dim] - curr_x[curr_dim];
        // dist2bb = dist2bb * dist2bb;
        dist2bbarr[curr_dim] = (query[curr_dim] - curr_x[curr_dim]) * (query[curr_dim] - curr_x[curr_dim]);
        T dist2bb = 0;
        for (int i = 0; i < dimension_; i++) {
            dist2bb += dist2bbarr[i];
        }
        if ((neighbor_idx.size() < k) || (dist2bb < distances.back())) {
            updateKNN<T>(query, curr_x, implicit_idx_tree_[node_idx], neighbor_idx, distances, k);
            searchKNNWithImplicitTreeOtherBranch(node_idx * 2 + 2, query, k, next_dim, neighbor_idx, distances, dist2bbarr);
        }
    } else {
        // advance to right in kdtree
        searchKNNWithImplicitTreeOtherBranch(node_idx * 2 + 2, query, k, next_dim, neighbor_idx, distances, dist2bbarr_cpy);
        // compute L2-distance^2 from the boundary to the query
        dist2bbarr[curr_dim] = (query[curr_dim] - curr_x[curr_dim]) * (query[curr_dim] - curr_x[curr_dim]);
        T dist2bb = 0;
        for (int i = 0; i < dimension_; i++) {
            dist2bb += dist2bbarr[i];
        }
        if ((neighbor_idx.size() < k) || (dist2bb < distances.back())) {
            updateKNN<T>(query, curr_x, implicit_idx_tree_[node_idx], neighbor_idx, distances, k);
            searchKNNWithImplicitTreeOtherBranch(node_idx * 2 + 1, query, k, next_dim, neighbor_idx, distances, dist2bbarr);
        }
    }

    free(dist2bbarr_cpy);

}

// a helper function for knn search
template<typename T>
void KDTree<T>::searchKNNWithImplicitTree(size_t node_idx,
                                          const std::vector<T>& query,
                                          const int k, int curr_dim,
                                          std::vector<size_t>& neighbor_idx,
                                          std::vector<T>& distances) {

    // index out of range of the implicit tree
    if (node_idx >= n_points_) {
        return;
    }

    // update number of nodes visited
    visited_++;

    // get current vector of interest
    const T* curr_x = data_ + implicit_idx_tree_[node_idx] * dimension_;

    int next_dim = (curr_dim + 1) % dimension_;

    if (query[curr_dim] <= curr_x[curr_dim]) {
        // advance to left in kdtree
        searchKNNWithImplicitTree(node_idx * 2 + 1, query, k, next_dim, neighbor_idx, distances);
        // compute L2-distance^2 from the boundary to the query
        T dist2bb = (query[curr_dim] - curr_x[curr_dim]) * (query[curr_dim] - curr_x[curr_dim]);
        if ((neighbor_idx.size() < k) || (dist2bb < distances.back())) {
            updateKNN<T>(query, curr_x, implicit_idx_tree_[node_idx], neighbor_idx, distances, k);
            T* dist2bbarr = (T*)malloc(sizeof(T)*dimension_);
            for (int i = 0; i < dimension_; i++) {
                dist2bbarr[i] = 0;
            }
            dist2bbarr[curr_dim] = dist2bb;
            searchKNNWithImplicitTreeOtherBranch(node_idx * 2 + 2, query, k, next_dim, neighbor_idx, distances, dist2bbarr);
        }
    } else {
        // advance to right in kdtree
        searchKNNWithImplicitTree(node_idx * 2 + 2, query, k, next_dim, neighbor_idx, distances);
        // compute L2-distance^2 from the boundary to the query
        T dist2bb = (query[curr_dim] - curr_x[curr_dim]) * (query[curr_dim] - curr_x[curr_dim]);
        if ((neighbor_idx.size() < k) || (dist2bb < distances.back())) {
            updateKNN<T>(query, curr_x, implicit_idx_tree_[node_idx], neighbor_idx, distances, k);
            T* dist2bbarr = (T*)malloc(sizeof(T)*dimension_);
            for (int i = 0; i < dimension_; i++) {
                dist2bbarr[i] = 0;
            }
            dist2bbarr[curr_dim] = dist2bb;
            searchKNNWithImplicitTreeOtherBranch(node_idx * 2 + 1, query, k, next_dim, neighbor_idx, distances, dist2bbarr);
        }
    }

}

// knn search
template<typename T>
int KDTree<T>::searchKNN(const std::vector<T>& query, const int k,
                         std::vector<size_t>& neighbor_idx,
                         std::vector<T>& distances) {

    // initialization
    visited_ = 0;
    neighbor_idx.clear();
    distances.clear();
    searchKNNWithImplicitTree(0, query, k, 0, neighbor_idx, distances);
    return visited_;
}

// a helper funciton for radius search
template <typename T>
void KDTree<T>::searchRadiusWithImplicitTree(size_t node_idx,
                                             const std::vector<T> &query,
                                             const T radius,
                                             int curr_dim,
                                             std::vector<size_t> &neighbor_idx,
                                             std::vector<T> &distances) {

    // index out of range of the implicit tree
    if (node_idx >= n_points_) {
        return;
    }

    // update number of nodes visited
    visited_++;
    // get current vector of interest
    const T* curr_x = data_ + implicit_idx_tree_[node_idx] * dimension_;

    // compute L2-distance^2 from the current vector to the query
    T distance = 0;
    for (int i = 0; i < dimension_; i++) {
        T a = curr_x[i] - query[i];
        distance += a * a;
    }

    // update k-neighbors
    if (distance <= radius) {
        distances.push_back(distance);
        neighbor_idx.push_back(implicit_idx_tree_[node_idx]);
        sortLastElementBasedOnDistance<T>(neighbor_idx, distances);
    }

    int next_dim = (curr_dim + 1) % dimension_;

    if (query[curr_dim] <= curr_x[curr_dim]) {
        // advance to left in kdtree
        searchRadiusWithImplicitTree(node_idx*2+1, query, radius, next_dim, neighbor_idx, distances);
        // compute L2-distance^2 from the boundary to the query
        T dist2bb = query[curr_dim] - curr_x[curr_dim];
        dist2bb = dist2bb * dist2bb;
        if (dist2bb < radius) {
            searchRadiusWithImplicitTree(node_idx*2+2, query, radius, next_dim, neighbor_idx, distances);
        }
    } else {
        // advance to right in kdtree
        searchRadiusWithImplicitTree(node_idx*2+2, query, radius, next_dim, neighbor_idx, distances);
        // compute L2-distance^2 from the boundary to the query
        T dist2bb = query[curr_dim] - curr_x[curr_dim];
        dist2bb = dist2bb * dist2bb;
        if (dist2bb < radius) {
            searchRadiusWithImplicitTree(node_idx*2+1, query, radius, next_dim, neighbor_idx, distances);
        }
    }

}

// radius search
template <typename T>
int KDTree<T>::searchRadius(const std::vector<T>& query, const T radius,
                            std::vector<size_t>& neighbor_idx,
                            std::vector<T>& distances) {

    // initialization
    visited_ = 0;
    neighbor_idx.clear();
    distances.clear();

    searchRadiusWithImplicitTree(0, query, radius, 0, neighbor_idx, distances);
    return visited_;
}

// a helper funciton for radius search
template <typename T>
void KDTree<T>::searchHybridWithImplicitTree(size_t node_idx,
                                             const std::vector<T> &query,
                                             const T radius,
                                             int curr_dim,
                                             std::vector<size_t> &neighbor_idx,
                                             std::vector<T> &distances,
                                             const size_t max_neighbors) {

    // index out of range of the implicit tree
    if (node_idx >= n_points_) {
        return;
    }

    // update number of nodes visited
    visited_++;
    // get current vector of interest
    const T* curr_x = data_ + implicit_idx_tree_[node_idx] * dimension_;

    // compute L2-distance^2 from the current vector to the query
    T distance = 0;
    for (int i = 0; i < dimension_; i++) {
        T a = curr_x[i] - query[i];
        distance += a * a;
    }

    // update k-neighbors
    if ((distance <= radius)&&(neighbor_idx.size() < max_neighbors)) {
        distances.push_back(distance);
        neighbor_idx.push_back(implicit_idx_tree_[node_idx]);
        sortLastElementBasedOnDistance<T>(neighbor_idx, distances);
    } else if (distance < distances.back()) {
        distances.back() = distance;
        neighbor_idx.back() = implicit_idx_tree_[node_idx];
        sortLastElementBasedOnDistance<T>(neighbor_idx, distances);
    }

    int next_dim = (curr_dim + 1) % dimension_;

    if (query[curr_dim] <= curr_x[curr_dim]) {
        // advance to left in kdtree
        searchHybridWithImplicitTree(node_idx*2+1, query, radius, next_dim, neighbor_idx, distances, max_neighbors);
        // compute L2-distance^2 from the boundary to the query
        T dist2bb = query[curr_dim] - curr_x[curr_dim];
        dist2bb = dist2bb * dist2bb;
        if (dist2bb < radius) {
            searchHybridWithImplicitTree(node_idx*2+2, query, radius, next_dim, neighbor_idx, distances, max_neighbors);
        }
    } else {
        // advance to right in kdtree
        searchHybridWithImplicitTree(node_idx*2+2, query, radius, next_dim, neighbor_idx, distances, max_neighbors);
        // compute L2-distance^2 from the boundary to the query
        T dist2bb = query[curr_dim] - curr_x[curr_dim];
        dist2bb = dist2bb * dist2bb;
        if (dist2bb < radius) {
            searchHybridWithImplicitTree(node_idx*2+1, query, radius, next_dim, neighbor_idx, distances, max_neighbors);
        }
    }

}

// radius search
template <typename T>
int KDTree<T>::searchHybrid(const std::vector<T>& query, const T radius,
                            std::vector<size_t>& neighbor_idx,
                            std::vector<T>& distances,
                            size_t max_neighbors /* = 0 */) {

    // initialization
    visited_ = 0;
    neighbor_idx.clear();
    distances.clear();
    if (max_neighbors == 0) {
        max_neighbors = n_points_;
    }

    searchHybridWithImplicitTree(0, query, radius, 0, neighbor_idx, distances, max_neighbors);
    return visited_;
}

template class KDTree<float>;
template class KDTree<double>;
