
#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "kdtree/kdtree.h"

// a helper function that performs insertion sort for the last element
template <typename T>
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

template <typename T>
void KDTree<T>::checkQueue(const std::vector<T>& query, const int k,
                           std::vector<size_t>& neighbor_idx,
                           std::vector<T>& distances) {

  if (queue_.empty()) {
     return;
  }

  NodeInfo next_node = queue_.top();
  queue_.pop();

  while (next_node.idx_ >= n_points_) {
      if (queue_.empty()) {
          return;
      }
      next_node = queue_.top();
      queue_.pop();
  }

  T dist2bb = next_node.distance_;

  if ((neighbor_idx.size() >= k)&&(distances.back() < dist2bb)) {
    return;
  }

  size_t node_idx = next_node.idx_;
  const T* curr_x = data_ + implicit_idx_tree_[node_idx] * dimension_;
  visited_++;
  updateKNN<T>(query, curr_x, implicit_idx_tree_[node_idx], neighbor_idx, distances, k);

  searchKNNWithImplicitTree(next_node.next_idx_, query, k, next_node.next_dim_,
                            neighbor_idx, distances, next_node.dist_vector_);

}

// a helper function for knn search
template <typename T>
void KDTree<T>::searchKNNWithImplicitTree(size_t node_idx,
                                          const std::vector<T>& query,
                                          const int k, int curr_dim,
                                          std::vector<size_t>& neighbor_idx,
                                          std::vector<T>& distances,
                                          std::vector<T>& dist2bbarr) {

    // index out of range of the implicit tree
    if (node_idx >= n_points_) {
        return;
    }

    // get current vector of interest
    const T* curr_x = data_ + implicit_idx_tree_[node_idx] * dimension_;

    std::vector<T> dist2bbarr_cpy = dist2bbarr;

    dist2bbarr[curr_dim] = (query[curr_dim] - curr_x[curr_dim]) * (query[curr_dim] - curr_x[curr_dim]);
    T dist2bb = 0;
    for (auto d : dist2bbarr) {
        dist2bb += d;
    }

    int next_dim = (curr_dim + 1) % dimension_;

    if (query[curr_dim] <= curr_x[curr_dim]) {
      queue_.emplace(node_idx, node_idx * 2 + 2, next_dim, dist2bb, dist2bbarr);
      searchKNNWithImplicitTree(node_idx * 2 + 1, query, k, next_dim, neighbor_idx, distances, dist2bbarr_cpy);
    } else {
      queue_.emplace(node_idx, node_idx * 2 + 1, next_dim, dist2bb, dist2bbarr);
      searchKNNWithImplicitTree(node_idx * 2 + 2, query, k, next_dim, neighbor_idx, distances, dist2bbarr_cpy);
    }

    //updateKNN<T>(query, curr_x, implicit_idx_tree_[node_idx], neighbor_idx, distances, k);
    checkQueue(query, k, neighbor_idx, distances);
}

// knn search
template<typename T>
int KDTree<T>::searchKNN(const std::vector<T>& query, const int k,
                         std::vector<size_t>& neighbor_idx,
                         std::vector<T>& distances) {

    if (k <= 0) {
       throw std::runtime_error("k must be bigger than 0.");
    }
    // initialization
    visited_ = 0;
    neighbor_idx.clear();
    distances.clear();
    std::vector<T> dist2bbarr(dimension_, 0);
    searchKNNWithImplicitTree(0, query, k, 0, neighbor_idx, distances, dist2bbarr);
    while (!queue_.empty()) {
      queue_.pop();
    }
    return visited_;
}

template <typename T>
static void updateRadius(const std::vector<T>& query,
                         const T* x, size_t x_idx,
                         std::vector<size_t>& neighbor_idx,
                         std::vector<T>& distances, double radius) {
    // compute L2-distance^2 from the current vector to the query
    T distance = 0;
    for (int i = 0; i < query.size(); i++) {
        T a = x[i] - query[i];
        distance += a * a;
    }

    // update k-neighbors
    if (distance <= radius) {
        distances.push_back(distance);
        neighbor_idx.push_back(x_idx);
        sortLastElementBasedOnDistance<T>(neighbor_idx, distances);
    }
}

// a helper funciton for radius search
template <typename T>
void KDTree<T>::searchRadiusWithImplicitTree(size_t node_idx,
                                             const std::vector<T>& query,
                                             const T radius,
                                             int curr_dim,
                                             std::vector<size_t>& neighbor_idx,
                                             std::vector<T>& distances,
                                             std::vector<T>& dist2bbarr) {

    // index out of range of the implicit tree
    if (node_idx >= n_points_) {
        return;
    }

    // get current vector of interest
    const T* curr_x = data_ + implicit_idx_tree_[node_idx] * dimension_;

    std::vector<T> dist2bbarr_cpy = dist2bbarr;

    dist2bbarr[curr_dim] = (query[curr_dim] - curr_x[curr_dim]) * (query[curr_dim] - curr_x[curr_dim]);
    T dist2bb = 0;
    for (auto d : dist2bbarr) {
        dist2bb += d;
    }

    int next_dim = (curr_dim + 1) % dimension_;

    if (query[curr_dim] <= curr_x[curr_dim]) {
        // advance to left in kdtree
        searchRadiusWithImplicitTree(node_idx*2+1, query, radius, next_dim, neighbor_idx, distances, dist2bbarr_cpy);
        // compute L2-distance^2 from the boundary to the query
        if (dist2bb <= radius) {
            visited_++;
            updateRadius<T>(query, curr_x, implicit_idx_tree_[node_idx], neighbor_idx, distances, radius);
            searchRadiusWithImplicitTree(node_idx*2+2, query, radius, next_dim, neighbor_idx, distances, dist2bbarr);
        }
    } else {
        // advance to right in kdtree
        searchRadiusWithImplicitTree(node_idx*2+2, query, radius, next_dim, neighbor_idx, distances, dist2bbarr_cpy);
        // compute L2-distance^2 from the boundary to the query
        if (dist2bb <= radius) {
            visited_++;
            updateRadius<T>(query, curr_x, implicit_idx_tree_[node_idx], neighbor_idx, distances, radius);
            searchRadiusWithImplicitTree(node_idx*2+1, query, radius, next_dim, neighbor_idx, distances, dist2bbarr);
        }
    }

}

// radius search
template <typename T>
int KDTree<T>::searchRadius(const std::vector<T>& query, const T radius,
                            std::vector<size_t>& neighbor_idx,
                            std::vector<T>& distances) {
    if (radius <= 0) {
        throw std::runtime_error("Radius must be bigger than 0.");
    }
    // initialization
    visited_ = 0;
    neighbor_idx.clear();
    distances.clear();
    std::vector<T> dist2bbarr(dimension_, 0);
    searchRadiusWithImplicitTree(0, query, radius * radius, 0, neighbor_idx, distances, dist2bbarr);
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
