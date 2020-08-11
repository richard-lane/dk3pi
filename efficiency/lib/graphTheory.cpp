#include <algorithm>

#include "graphTheory.h"

Edge::Edge(const size_t node, const double weight) : _weight(weight), _node(node)
{
    ;
}

bool Edge::operator==(const Edge& other) const
{
    return _weight == other._weight;
}

bool Edge::operator!=(const Edge& other) const
{
    return !(*this == other);
}

bool Edge::operator<(const Edge& other) const
{
    return _weight < other._weight;
}

bool Edge::operator>(const Edge& other) const
{
    return !(*this < other);
}

Graph::Graph(const size_t numNodes) : _numNodes(numNodes), _adjacencyList(std::vector<std::list<Edge>>(_numNodes))
{
    ;
}

void Graph::addEdge(const size_t node1, const size_t node2, const double weight)
{
    _adjacencyList[node1].push_back(Edge(node2, weight));
    _adjacencyList[node2].push_back(Edge(node1, weight));
}

bool Graph::_makesCycle(const size_t node1, const size_t node2) const
{
    // Copy our adjacency list
    std::vector<std::list<Edge>> adjList = _adjacencyList;

    // Add our edge to the new adjacency list with an arbitrary weight of 1
    // This is repeated code, but i don't care
    adjList[node1].push_back(Edge(node2, 1));
    adjList[node2].push_back(Edge(node1, 1));

    // Perform a Depth first search on the graph for any cycles
    // Start the search at the one of the new nodes, since we know that if a cycle has been newly formed it must contain
    // this node
    std::vector<bool> discovered(_numNodes, false);
    return containsCycle(adjList, node1, discovered, node2);
}

bool containsCycle(const std::vector<std::list<Edge>>& adjacencyList,
                   const size_t                        node,
                   std::vector<bool>&                  discovered,
                   const size_t                        parent)
{
    // Could do a short circuit-test here; if there are more edges than nodes, then there is definitely a cycle

    // Mark this node as discovered
    discovered[node] = true;

    // Look at each node connected to the one passed in to this function
    for (Edge edge : adjacencyList[node]) {

        // If this node has not been discovered yet, start the depth-first search again on this node
        if (!discovered[edge.node()]) {
            if (containsCycle(adjacencyList, edge.node(), discovered, node)) {
                return true;
            }

            // If this node has been discovered and is not our node's parent, we have a cycle
        } else if (edge.node() != parent) {
            return true;
        }
    }

    // No back-edges found
    return false;
}

bool isSpanning(const std::vector<std::list<Edge>>& adjacencyList)
{
    // Could do a short-circuit test here: if there are fewer edges than (nodes-1), then our graph cannot be spanning

    // Repeated code here from the containsCycle function, but i don't care again
    // Start at node 0
    std::vector<bool> discoveredNodes(adjacencyList.size(), false);
    bool              newNodeFound{true};
    discoveredNodes[0] = true;
    while (newNodeFound) {
        newNodeFound = false;
        // Iterate over our vector of discovered nodes
        for (size_t i = 0; i < discoveredNodes.size(); ++i) {
            // If a node is discovered, find all the nodes it connects to
            if (discoveredNodes[i]) {
                for (Edge edge : adjacencyList[i]) {
                    if (!discoveredNodes[edge.node()]) {
                        newNodeFound                 = true;
                        discoveredNodes[edge.node()] = true;
                    }
                }
            }
        }
    }
    return std::all_of(discoveredNodes.begin(), discoveredNodes.end(), [](const bool x) { return x; });
}
