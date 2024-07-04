#ifndef COMMON_H
#define COMMON_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <Eigen/Dense>
#include <vector>
#include <unordered_set>
#include <algorithm> 

namespace py = pybind11;

using MatrixXs = Eigen::Matrix<unsigned short, Eigen::Dynamic, Eigen::Dynamic>;
using VectorXs =Eigen::Matrix<unsigned short, Eigen::Dynamic, 1>;

using Vectors = Eigen::Map<const Eigen::Matrix<unsigned short, Eigen::Dynamic, 1>>;
using Dict_cpp = std::unordered_map<int, MatrixXs>;
using Array_py =  py::array_t<unsigned short, 2>;
using Dict_py = std::unordered_map<int, py::array_t<unsigned short, 2>>;



#endif // COMMON_H