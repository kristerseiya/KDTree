
#pragma once

#include <cstdlib>
#include <vector>
#include <cstdbool>

template <typename T>
class KDTree {

protected:
    size_t n_points_;
    int dimension_;
    size_t visited_;
    const T* data_;
    size_t* implicit_idx_tree_;
    bool copied_;

public:
  KDTree();
  KDTree(const T* data, int dim, size_t num, bool copy = false);
  KDTree(const std::vector<T>& data, int dim, bool copy = false);
  ~KDTree();

  void assign(const T* data, int dim, size_t num, bool copy = false);
  void assign(const std::vector<T>& data, int dim, bool copy = false);
  void add(const T* data, size_t num);

  bool isCopied();

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
  void searchKNNWithImplicitTreeOtherBranch(size_t node_idx,
                                 const std::vector<T>& query,
                                 const int k, int curr_dim,
                                 std::vector<size_t>& neighbor_idx,
                                 std::vector<T>& distances,
                                 T* dist2bbarr);
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
