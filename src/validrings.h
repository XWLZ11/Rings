#ifndef VALIDRINGS_H
#define VALIDRINGS_H

#include "py2cpp.h"

bool new_bond(py::array_t<unsigned short> ring_py, const Dict_py& cycles_single_py, const Dict_py& ref_rings_py, const int& ring_definition);
std::vector<bool> ring_in_rings(const VectorXs& ring, std::vector<bool>& bond_exists, const MatrixXs& rings);
bool shortcutted_by_ring(const VectorXs& ring, const MatrixXs& rings);
unsigned short count_overlap(std::vector<unsigned short> ring, const VectorXs& tmp_vec);



class RingsGroup {
public:
    RingsGroup(const Dict_cpp& rings);
    ~RingsGroup();
	Dict_cpp validing(const int& ring_definition);

private:
	Dict_cpp rings_all;
	Dict_cpp rings_single;
	Dict_cpp rings_ref;
	
	template<typename T>
	bool new_bond(const T& ring, const int& ring_definition);
};






#endif   // VALIDRINGS_H