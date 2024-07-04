#include "findrings.h"
#include "validrings.h"
#include <iostream>


py::list findrings(const py::array_t<unsigned short> bonds_py, const unsigned short& num_nodes, const unsigned short& length_bound) {
    Eigen::Map<MatrixXs> bonds = array_py2cpp(bonds_py);
    Graph graph(bonds, num_nodes);    
    std::vector<Ring> rings = graph.findAllRings(length_bound);
    return rings_to_list(rings);
}

Dict_py validings(const Dict_py& sum_rings_py, const int& ring_definition) {
    py::gil_scoped_acquire acquire; // 获取 GIL

    Dict_cpp sum_rings_cpp = dict_py2cpp(sum_rings_py);
    RingsGroup rings_group(sum_rings_cpp);
    Dict_cpp rings_single_cpp = rings_group.validing(ring_definition);
    Dict_py rings_single_py = dict_cpp2py(rings_single_cpp);
    
    /**
     * 找了很长时间的bug，发现这里传递数据存在问题，必须输出一下才可以正常运行，
     * 不然python会报类型错误，目前不知道std::cout对数据输出前进行了怎样的初
     * 始化操作，所以暂时先禁止输出，后续有时间再研究一下。
     * 可能存在的问题：
     * 1. python的GIL问题，目前我只会获取和释放，比较难排查
     * 2. 输出前的初始化问题
     * 3. 并发执行问题，可能性较大
     * 3. 类型注册问题，可能需要在python中注册一下类型，但我感觉可能小比较小
     */
    
    std::cout.setstate(std::ios::failbit); // 禁用输出
    std::cout << rings_single_py[0] << std::endl;
    std::cout.clear(); // 恢复输出

    py::gil_scoped_release release; // 释放 GIL
    return rings_single_py;
}


PYBIND11_MODULE(rings, m) {
    m.def("findrings", &findrings);
    m.def("validrings", &validings);   
    m.def("new_bond", &new_bond);
}