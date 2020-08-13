#include <algorithm>
#include <functional>

#include "graphTheory.h"

Graph::Graph(const size_t numNodes) : _numNodes(numNodes), _adjacencyList(std::vector<std::list<Edge>>(_numNodes))
{
    ;
}

void Graph::addEdge(const size_t node1, const size_t node2, const double weight)
{
    // Should possibly also check that the edge doesn't already exist

    _adjacencyList[node1].push_back(Edge(node1, node2, weight));
    _adjacencyList[node2].push_back(Edge(node2, node1, weight));
}

std::vector<std::list<Edge>> Graph::getMaxSpanningTree() const
{
    // Create a list of our edges in order of weight
    std::vector<Edge> edges;
    size_t            numEdges{0};
    for (auto list : _adjacencyList) {
        numEdges += list.size();
    }
    edges.reserve(numEdges);
    for (auto list : _adjacencyList) {
        edges.insert(edges.end(), list.begin(), list.end());
    }
    std::sort(edges.begin(), edges.end(), std::greater<Edge>());

    std::vector<std::list<Edge>> mst(_numNodes);
    for (Edge edge : edges) {
        if (!makesCycle(edge.from(), edge.to(), mst)) {
            mst[edge.from()].push_back(edge);
        }
        if (isSpanning(mst)) {
            break; // don't like it
        }
    }
    return mst;
}

bool makesCycle(const size_t node1, const size_t node2, const std::vector<std::list<Edge>>& adjacencyList)
{
    // Copy our adjacency list
    std::vector<std::list<Edge>> adjList = adjacencyList;

    // Add our edge to the new adjacency list with an arbitrary weight of 1
    // This is repeated code, but i don't care
    adjList[node1].push_back(Edge(node1, node2, 1));
    adjList[node2].push_back(Edge(node2, node1, 1));

    // Perform a Depth first search on the graph for any cycles
    // Start the search at the one of the new nodes, since we know that if a cycle has been newly formed it must contain
    // this node
    const size_t      numNodes{adjList.size()};
    std::vector<bool> discovered(numNodes, false);
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
        if (!discovered[edge.to()]) {
            if (containsCycle(adjacencyList, edge.to(), discovered, node)) {
                return true;
            }

            // If this node has been discovered and is not our node's parent, we have a cycle
        } else if (edge.to() != parent) {
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
                    if (!discoveredNodes[edge.to()]) {
                        newNodeFound               = true;
                        discoveredNodes[edge.to()] = true;
                    }
                }
            }
        }
    }
    return std::all_of(discoveredNodes.begin(), discoveredNodes.end(), [](const bool x) { return x; });
}

std::vector<std::list<Edge>> inTree(const size_t root, const std::vector<std::list<Edge>>& treeAdjacencyList)
{
    size_t            order{treeAdjacencyList.size()};
    std::vector<bool> discovered(order, false);
    if (containsCycle(treeAdjacencyList, root, discovered)) {
        throw CyclicGraph();
    }

    if (!isSpanning(treeAdjacencyList)) {
        throw NotSpanning();
    }

    // Reset our vector of discovered nodes to false, except the one we start at
    std::fill(discovered.begin(), discovered.end(), false);
    discovered[root] = true;

    std::vector<std::list<Edge>> directedTree(order);

    // Iterate until we've found all our vertices
    while (!std::all_of(discovered.begin(), discovered.end(), [](bool v) { return v; })) {

        // Iterate over our vertices
        for (size_t i = 0; i < order; ++i) {
            for (Edge edge : treeAdjacencyList[i]) {
                if (!discovered[i] && discovered[edge.to()]) {
                    // Draw an edge from this vertex to the discovered vertex
                    directedTree[i].push_back(edge);
                    discovered[edge.from()] = true;
                }
            }
        }
    }

    return directedTree;
}
