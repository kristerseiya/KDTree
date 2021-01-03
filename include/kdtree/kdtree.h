
#pragma once

#include <cstdlib>
#include <vector>

template <typename T>
class KDTree {

public:
    size_t num_nodes_;
    int dimension_;
    size_t visited_;
protected:
    const T* data_;
    size_t* implicit_idx_tree_;

public:
  KDTree();
  KDTree(const T* data, int dim, size_t num);
  KDTree(const std::vector<T>& data, int dim);
  ~KDTree();

  void set(const T* data, int dim, size_t num);
  int searchKNN(const std::vector<T>& query, const int k,
                std::vector<size_t>& neighbor_idx,
                std::vector<T>& distances);
  int searchRadius(const std::vector<T>& query, const T radius,
                   std::vector<size_t>& neighbor_idx,
                   std::vector<T>& distances);
  int searchHybrid(const std::vector<T>& query, const T radius,
                     std::vector<size_t>& neighbor_idx,
                     std::vector<T>& distances,
                     size_t max_neighbors = 0);

private:
  void searchKNNWithImplicitTree(size_t node_idx,
                                 const std::vector<T>& query,
                                 const int k, int curr_dim,
                                 std::vector<size_t>& neighbor_idx,
                                 std::vector<T>& distances);
  void searchRadiusWithImplicitTree(size_t node_idx,
                                    const std::vector<T> &query,
                                    const T radius,
                                    int curr_dim,
                                    std::vector<size_t> &neighbor_idx,
                                    std::vector<T> &distances);
  void searchHybridWithImplicitTree(size_t node_idx,
                                    const std::vector<T> &query,
                                    const T radius,
                                    int curr_dim,
                                    std::vector<size_t> &neighbor_idx,
                                    std::vector<T> &distances,
                                    const size_t max_neighbors);
};

