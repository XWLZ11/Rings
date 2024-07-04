#include "findrings.h"

#include <iostream>

#include "py2cpp.h"




Ring::Ring(const std::vector<unsigned short>& nodes) : nodes(nodes) {
}
Ring::~Ring() {
}
bool Ring::operator==(const Ring& other) const {
    return hash == other.hash;
}
size_t Ring::getHash() const {
    std::vector<unsigned short> sorted_nodes = nodes; 
    sort(sorted_nodes.begin(), sorted_nodes.end());
    size_t seed = 0;
    for (unsigned short node : sorted_nodes) { 
        seed ^= std::hash<unsigned short>{}(node) + 0x9e3779b9 + (seed << 6) + (seed >> 2);   
    } 
    return seed;
}

Graph::Graph(const Eigen::Map<MatrixXs>& edges, const unsigned short& num_nodes) : neighbors(num_nodes) {
    for (int i = 0; i < edges.rows(); ++i) {
        this -> addEdge(edges(i, 0), edges(i, 1));
    }
}

Graph::~Graph() {
}
void Graph::addEdge(unsigned short u, unsigned short v) {
    neighbors[u].push_back(v);
    neighbors[v].push_back(u); // 无向图需要双向连接
}


std::vector<Ring> Graph::findAllRings(unsigned short size_limit) {
    
    std::vector<Ring> rings;
    unsigned short num_nodes = neighbors.size(); 
    std::vector<bool> visited(num_nodes, false); 
    // 遍历图中的每个节点，从每个节点开始深度优先搜索
    for (unsigned short i = 0; i < num_nodes; ++i) {
        std::vector<unsigned short> path;
        dfs(i, i, size_limit, path, visited, rings);
    }
    return rings;
}

void Graph::dfs(unsigned short current, unsigned short start, unsigned short size_limit, std::vector<unsigned short>& path, std::vector<bool>& visited, std::vector<Ring>& rings) {
    visited.at(current) = true; // 标记当前节点为已访问
    path.push_back(current); 
    unsigned short current_size = path.size();
    for (unsigned short neighbor : neighbors[current]) {
        if (neighbor == start && current_size >= 3 && current_size <= size_limit) { // 找到一个环
            Ring ring(path);
            if (rings_set.find(ring) == rings_set.end()) {
                rings.push_back(ring);
                rings_set.insert(ring);
            }
            continue;
        }            
        if (!visited[neighbor] && current_size <= size_limit) { 
            dfs(neighbor, start, size_limit, path, visited, rings);
        }            
    }
    visited.at(current) = false; // 回溯前将节点标记为未访问
    path.pop_back(); // 回溯前移除当前节点
}

