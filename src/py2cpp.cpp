#include "py2cpp.h"


std::vector<unsigned short> array2vector(const py::array_t<unsigned short>& array) {
    size_t* ptr = static_cast<size_t*>(array.request().ptr);
    std::vector<unsigned short> vec(ptr, ptr + array.size());
    return vec;
}


Eigen::Map<MatrixXs> array_py2cpp(const Array_py& array_py) {
    // 这里读取进来的数据是以列优先的，即先读一列再读下一列，与python传入的array行存储不同，需要注意
    py::buffer_info buf = array_py.request(); 
    // 这里的类型要与python传入的bonds数据类型一致（unsigned short == np.uint16）
    unsigned short* ptr = static_cast<unsigned short*>(buf.ptr); 
    
    size_t nrows = buf.shape[0];
    unsigned short ncols = buf.shape[1];
    
    Eigen::Map<MatrixXs> matrix_cpp(ptr, nrows, ncols);
    /** 用数组速度可能会快一点（与vector不会相差很大），但会使代码变得臃肿，
     * 因为要在调用前创建数组，传入也得用指针或者引用，不如直接用vector更方便，
     * 用Eigen的矩阵是因为可以直接映射而不必赋值，总之，这里对于效率影响微乎其微
    **/
    return matrix_cpp;
}
Array_py array_cpp2py(const MatrixXs& array_cpp) {
    Array_py array_py = py::array_t<unsigned short>(
        {array_cpp.rows(), array_cpp.cols()}, 
        {sizeof(unsigned short) * array_cpp.cols(), 
         sizeof(unsigned short)}, 
         array_cpp.data(),
         py::cast(array_cpp)
    ); 
    return array_py;
}

Dict_cpp dict_py2cpp(const Dict_py& dict_py) {
	Dict_cpp dict_cpp;
    // 遍历每个key-value对
	for (const auto& entry : dict_py) { 
        dict_cpp[entry.first] = array_py2cpp(entry.second).transpose();
	}
	return dict_cpp;
}

Dict_py dict_cpp2py(const Dict_cpp& dict_cpp){
    Dict_py dict_py;
    for (const auto& entry : dict_cpp) {
        dict_py[entry.first] = array_cpp2py(entry.second.transpose());
    }
    return dict_py;
}

py::list rings_to_list(const std::vector<Ring>& rings) {
    py::list lists;
    for (const Ring& ring : rings) {
        py::list nodes_list;
        for (unsigned short node : ring.nodes) {
            nodes_list.append(node);
        }
        lists.append(nodes_list);
    }
    return lists;
}



	