//
// Created by Krister Ulvog on 1/1/21.
//

#pragma once

#include <Eigen/StdVector>
#include <string>

void ReadPLY(const std::string& filename, std::vector<Eigen::Vector3d>& points,
             std::vector<Eigen::Vector3d>& normals,
             std::vector<Eigen::Vector3d>& colors);
