#ifndef FINDRINGS_H
#define FINDRINGS_H

#include "common.h"


class Ring {
public:
    Ring(const std::vector<unsigned short>& nodes);
    ~Ring();

    std::vector<unsigned short> nodes;
    size_t hash = getHash();
    
    bool operator==(const Ring& other) const;

private:
    size_t getHash() const;
};

// 自定义哈希函数

namespace std {
    template <>
    struct hash<Ring> {
        size_t operator() (const Ring& ring) const {   
            return ring.hash;
        }
    };
}


class Graph {
public:
    Graph(const Eigen::Map<MatrixXs>& edges, const unsigned short& num_nodes);
    ~Graph();
    void addEdge(unsigned short u, unsigned short v);
    std::vector<Ring> findAllRings(unsigned short size_limit);

private:
	std::vector<std::vector<unsigned short>> neighbors;
    std::unordered_set<Ring> rings_set; 
    
    void dfs(unsigned short current, unsigned short start, unsigned short size_limit, std::vector<unsigned short>& path, std::vector<bool>& visited, std::vector<Ring>& rings);
};

#endif   // FINDRINGS_H