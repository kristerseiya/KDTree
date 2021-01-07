#include <kdtree/kdtree.h>
#include "load_ply.h"
#include <Eigen/Core>
#include <iostream>

int main(int argc, char** argv) {
  if (argc < 2) {
    throw std::runtime_error("Not Enough Arguments");
  }

  std::vector<Eigen::Vector3d> points;
  std::vector<Eigen::Vector3d> normal;
  std::vector<Eigen::Vector3d> color;

  ReadPLY(argv[1], points,
          normal,
          color);

  std::cout << points.size() << std::endl;

  auto point_arr = Eigen::Map<const Eigen::MatrixXd>((const double*)points.data(), 1, points.size() * 3);
  KDTree<double> kdt(point_arr.data(),3,points.size(), false);

  std::vector<double> query(3);
  query[0] = points[0][0]; query[1] = points[0][1]; query[2] = points[0][2];
  std::cout << "My Query" << std::endl;
  std::cout << query[0] << ", " << query[1] << ", " << query[2] << std::endl;

  std::vector<size_t> neighbor_idx;
  std::vector<double> distances;
  kdt.searchKNN(query,
                20,
                neighbor_idx,
                distances);

  std::cout << std::endl;
  std::cout << "20 nearest neighbors" << std::endl;
  for (int i = 0; i < 20; i++) {
      std::cout << points[neighbor_idx[i]].transpose() << ": " << distances[i] << std::endl;
  }

  double radius = distances[9];
  kdt.searchRadius(query,
                   radius,
                   neighbor_idx,
                   distances);

  std::cout << std::endl;
  std::cout << radius << " radius search" << std::endl;
    for (int i = 0; i < neighbor_idx.size(); i++) {
        std::cout << points[neighbor_idx[i]].transpose() << ": " << distances[i] << std::endl;
    }
  return 0;
}
