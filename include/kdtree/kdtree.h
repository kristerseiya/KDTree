
#pragma once

#include <cstdlib>
#include <vector>
#include <cstdbool>
// #include <tuple>
#include <queue>

template <typename T>
class KDTree {

protected:
    size_t n_points_;
    int dimension_;
    size_t visited_;
    const T* data_;
    size_t* implicit_idx_tree_;
    bool copied_;

private:

    // struct NodeInfo {
    //   size_t idx_;
    //   int dim_;
    //   T distance_;
    //   std::vector<T> dist_vector_;
    //   NodeInfo(size_t idx, int dim, T distance, std::vector<T>& dist_vector) :
    //       idx_(idx), dim_(dim), distance_(distance), dist_vector_(dist_vector) {}
    // };
    struct NodeInfo {
      size_t idx_;
      size_t next_idx_;
      int next_dim_;
      T distance_;
      std::vector<T> dist_vector_;
      NodeInfo(size_t idx, size_t next_idx, int next_dim, T distance, std::vector<T>& dist_vector) :
          idx_(idx), next_idx_(next_idx), next_dim_(next_dim), distance_(distance), dist_vector_(dist_vector) {}
    };

    struct CompareNodeInfo {
      bool operator()(const NodeInfo& lhs, const NodeInfo& rhs) {
        return lhs.distance_ > rhs.distance_;
      }
    };

    std::priority_queue<NodeInfo, std::vector<NodeInfo>, CompareNodeInfo> queue_;

public:
  KDTree();
  KDTree(const T* data, int dim, size_t num, bool copy = false);
  KDTree(const std::vector<T>& data, int dim, bool copy = false);
  ~KDTree();

  void assign(const T* data, int dim, size_t num, bool copy = false);
  void assign(const std::vector<T>& data, int dim, bool copy = false);
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

  void checkQueue(const std::vector<T>& query, const int k,
                  std::vector<size_t>& neighbor_idx,
                  std::vector<T>& distances);

  void searchKNNWithImplicitTree(size_t node_idx,
                                 const std::vector<T>& query,
                                 const int k, int curr_dim,
                                 std::vector<size_t>& neighbor_idx,
                                 std::vector<T>& distances,
                                 std::vector<T>& dist2bbarr);
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
