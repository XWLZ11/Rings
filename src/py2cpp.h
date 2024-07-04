#ifndef PY2CPP_H
#define PY2CPP_H

#include "common.h"
#include "findrings.h"

std::vector<unsigned short> array2vector(const py::array_t<unsigned short>& array);

Eigen::Map<MatrixXs> array_py2cpp(const Array_py& array_py);
Array_py array_cpp2py(const MatrixXs& array_cpp);

py::list rings_to_list(const std::vector<Ring>& rings);

Dict_cpp dict_py2cpp(const Dict_py& dict_py);
Dict_py dict_cpp2py(const Dict_cpp& dict_cpp);


#endif   // PY2CPP_H