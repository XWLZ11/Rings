#include "validrings.h"




bool new_bond(py::array_t<unsigned short> ring_py, const Dict_py& cycles_single_py, const Dict_py& ref_rings_py, const int& ring_definition) {

	bool flag = true;

	Vectors ring(ring_py.data(), ring_py.size());

	if (ring_definition == 1) {
		unsigned short ring_size = ring.size();
		Dict_cpp cycles_single = dict_py2cpp(cycles_single_py);
		std::vector<bool> bond_exists(ring_size, false);
		for (const auto& entry : cycles_single) {		
			const MatrixXs& rings = entry.second;
			bond_exists = ring_in_rings(ring, bond_exists, rings);
		}
		flag = std::any_of(bond_exists.begin(), bond_exists.end(), [](bool b) { return !b; });
	}
	else if (ring_definition == 2) {
		Dict_cpp ref_rings = dict_py2cpp(ref_rings_py);
		for (const auto& entry : ref_rings) {
			const MatrixXs& rings = entry.second;
			if (shortcutted_by_ring(ring, rings)) {
				flag = false;
				break;
			}
		}
	}
	return flag;
}



std::vector<bool> ring_in_rings(const VectorXs& ring, std::vector<bool>& bond_exists, const MatrixXs& rings) {
	unsigned short ring_size = ring.size(), nrows = rings.rows(), ncols = rings.cols();	

	for (Eigen::Index col = 0; col < ncols; ++col) {
		for (unsigned short ia = 0; ia < ring_size; ++ia) {
			unsigned short i1 = ia, i2 = (ia + 1) % ring_size, i3 = (ia + 2) % ring_size;
			for (unsigned short ib = 0; ib < nrows; ++ib) {
				unsigned short j1 = ib, j2 = (ib + 1) % nrows, j3 = (ib + 2) % nrows;
				if ((rings(j1, col) == ring(i1) && 
					rings(j2, col) == ring(i2) && 
					rings(j3, col) == ring(i3)) ||
				   (rings(j1, col) == ring(i3) && 
					rings(j2, col) == ring(i2) &&
					rings(j3, col) == ring(i1))) {

					bond_exists.at(ia) = true;
					break;
				}
			}
		}
	} 
	return bond_exists;
}



bool is_category_1(const unsigned short& na_ref_ring, const unsigned short& nover) {
    return (na_ref_ring == 3 && nover == 3);
}
bool is_category_2(const unsigned short& na_ref_ring, const unsigned short& nover) {
    return ((na_ref_ring == 4 && nover == 4) ||
		    (na_ref_ring == 5 && nover == 4));
}
bool is_category_3(const unsigned short& na_ref_ring, const unsigned short& nover) {
    return ((na_ref_ring == 5 && nover == 5) ||
			(na_ref_ring == 6 && nover == 5) ||
		    (na_ref_ring == 7 && nover == 5));
}
bool is_category_4(const unsigned short& na_ref_ring, const unsigned short& nover) {
    return ((na_ref_ring == 6 && nover == 6) ||
			(na_ref_ring == 7 && nover == 6) ||
			(na_ref_ring == 8 && nover == 6) ||
		    (na_ref_ring == 9 && nover == 6));
}
bool is_category_5(const unsigned short& na_ref_ring, const unsigned short& nover) {
    return ((na_ref_ring == 7 && nover == 7) ||
			(na_ref_ring == 8 && nover == 7) ||
			(na_ref_ring == 9 && nover == 7) ||
			(na_ref_ring == 10 && nover == 7) ||
		    (na_ref_ring == 11 && nover == 7));
}
bool is_category_6(const unsigned short& na_ref_ring, const unsigned short& nover) {
    return ((na_ref_ring == 8 && nover == 8) ||
			(na_ref_ring == 9 && nover == 8) ||
			(na_ref_ring == 10 && nover == 8) ||
			(na_ref_ring == 11 && nover == 8) ||
			(na_ref_ring == 12 && nover == 8) ||
		    (na_ref_ring == 13 && nover == 8));
}
bool is_category_7(const unsigned short& na_ref_ring, const unsigned short& nover) {
    return ((na_ref_ring == 9 && nover == 9) ||
			(na_ref_ring == 10 && nover == 9) ||
			(na_ref_ring == 11 && nover == 9) ||
			(na_ref_ring == 12 && nover == 9) ||
			(na_ref_ring == 13 && nover == 9) ||
			(na_ref_ring == 14 && nover == 9) ||
		    (na_ref_ring == 15 && nover == 9));
}
bool is_category_8(const unsigned short& na_ref_ring, const unsigned short& nover) {
    return ((na_ref_ring == 10 && nover == 10) ||
			(na_ref_ring == 11 && nover == 10) ||
			(na_ref_ring == 12 && nover == 10) ||
			(na_ref_ring == 13 && nover == 10) ||
			(na_ref_ring == 14 && nover == 10) ||
			(na_ref_ring == 15 && nover == 10) ||
			(na_ref_ring == 16 && nover == 10) ||
		    (na_ref_ring == 17 && nover == 10));
}
bool shortcutted_by_ring(const VectorXs& ring, const MatrixXs& rings) {
	unsigned short ring_size = ring.size(), nrows = rings.rows(), ncols = rings.cols();
	bool flag = false;
	std::vector<unsigned short> cycle_vec(ring.data(), ring.data() + ring.size());
	for (Eigen::Index col = 0; col < ncols; ++col) {
		Vectors ref_ring_tmp(rings.col(col).data(), nrows);
		unsigned short nover = count_overlap(cycle_vec, ref_ring_tmp);
		switch (ring_size) {
			case 19:
			case 18:
				flag = flag || is_category_8(nrows, nover);
			case 17:
			case 16:
				flag = flag || is_category_7(nrows, nover);
			case 15:
			case 14:
				flag = flag || is_category_6(nrows, nover);
			case 13:
			case 12:
				flag = flag || is_category_5(nrows, nover);
			case 11:
			case 10:
				flag = flag || is_category_4(nrows, nover);
			case 9:
			case 8:
				flag = flag || is_category_3(nrows, nover);
			case 7:
			case 6:
				flag = flag || is_category_2(nrows, nover);
			case 5:
			case 4:
				flag = flag || is_category_1(nrows, nover);	
				break;

            default:
                continue;
        }
		if (flag) {
			break;
		}
	}
	return flag;
}

unsigned short count_overlap(std::vector<unsigned short> ring, const VectorXs& tmp_vec) {
	std::vector<unsigned short> ref_ring(tmp_vec.data(), tmp_vec.data() + tmp_vec.size());
	std::sort(ring.begin(), ring.end());
	std::sort(ref_ring.begin(), ref_ring.end());
	std::vector<unsigned short> intersection;
	std::set_intersection(ring.begin(), ring.end(),
                          ref_ring.begin(), ref_ring.end(),
                          std::back_inserter(intersection));
	return intersection.size();
}









template<typename T>
bool RingsGroup::new_bond(const T& ring, const int& ring_definition) {

	bool flag = true;

	if (ring_definition == 1) {
		unsigned short ring_size = ring.size();
		std::vector<bool> bond_exists(ring_size, false);
		for (const auto& entry : this -> rings_single) {		
			const MatrixXs& rings = entry.second;
			bond_exists = ring_in_rings(ring, bond_exists, rings);
		}
		flag = std::any_of(bond_exists.begin(), bond_exists.end(), [](bool b) { return !b; });
	}
	else if (ring_definition == 2) {
		for (const auto& entry : this -> rings_ref) {
			const MatrixXs& rings = entry.second;
			if (shortcutted_by_ring(ring, rings)) {
				flag = false;
				break;
			}
		}
	}
	return flag;
}



RingsGroup::RingsGroup(const Dict_cpp& rings) : rings_all(rings) {
	MatrixXs empty_matrix;
	size_t key_count = this -> rings_all.size();
	for (size_t i = 3; i < 3 + key_count; ++i) {
		this -> rings_single[i] = empty_matrix;
		this -> rings_ref[i] = empty_matrix;
	}
}
RingsGroup::~RingsGroup() {
}

Dict_cpp RingsGroup::validing(const int& ring_definition) {
	
	std::vector<VectorXs> empty_vector, tmp;
	for (const auto& entry : this -> rings_all) {
		tmp = empty_vector;
		for (Eigen::Index i = 0; i < entry.second.cols(); ++i) {
			const VectorXs& ring =  entry.second.col(i);
			if (this -> new_bond(ring, ring_definition)) {
				tmp.push_back(ring);
			}
		}
		this -> rings_ref[entry.first] = entry.second;
		if (tmp.size() == 0) {
			continue;
		}
		MatrixXs tmp_matrix(entry.first, tmp.size());
		for (size_t i = 0; i < tmp.size(); ++i) {
			tmp_matrix.col(i) = tmp[i];
		}
		this -> rings_single[entry.first] = tmp_matrix;
	}
	return this -> rings_single;
}



