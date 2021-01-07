# KDTree
KDTree for K Nearest Search and Radius Search in C++ and Python. **Faster build time than scipy.spatial.KDTree!**

## Building and Installing
Use CMAKE to build the library. \
First, edit CMakeLists.txt for configuration.
```bash
mkdir build && cd build && cmake ..
```
No 3rd party libraries required for C++ \
Pybind11 is required for building python interface.

## Usage in C++
See include/kdtree/kdtree.h to see the structure of KDTree class. \
Construct KDTree object with or without the data array
```C++
KDTree<double>::KDTree();
KDTree<double>::KDTree(const double* data, int dim, size_t num);
KDTree<double>::KDTree(const std::vector<double>& data, int dim);
```
To change or assign data after construction
```C++
KDTree<double>::assign(const double* data, int dim, size_t num);
```
You can also use ```KDTree<float> ``` if you want to use ```float``` type instead. \
To perform search, run one of the following member functions. \
These functions will return the number of nodes visited during the search. \
**Pass only one query at a time.**
```C++
// K Nearest Neighbor Search
int KDTree<double>::searchKNN(const std::vector<double>& query,
                              const int k,
                              std::vector<size_t>& neighbor_idx,
                              std::vector<double>& distances);
// Radius Search
int KDTree<double>::searchRadius(const std::vector<double>& query,
                                 const double radius,
                                 std::vector<size_t>& neighbor_idx,
                                 std::vector<double>& distances);
// Both K Nearest and Radius Search
int KDTree<double>::searchHybrid(const std::vector<double>& query,
                                 const double radius,
                                 std::vector<size_t>& neighbor_idx,
                                 std::vector<double>& distances,
                                 size_t max_neighbors = 0);
```
```std::vector<size_t> neighbor_idx``` will store the indices of the data. \
```std::vector<double> distances``` will store the distances from the query.
``` C++
std::vector<size_t> neighbor_idx // indices in the data (idx * dimension)
std::vector<double> distances // distances from the query
```
## Usage in Python
```Python
from pykdtree import KDTree

kdt = KDTree()
kdt.assign('''n x m numpy array''')
idx, dist = kdt.searchKNN(query, k :int)
idx, dist = kdt.searchRadius(query, radius :double)
idx, dist = kdt.searchHyBrid(query, radius :double, max_k :int)
```
