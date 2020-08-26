#include <algorithm>
#include <functional>
#include <stack>

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
    // edges will be an array of edges in order of weight- these will be e.g. {(0, 1), (1, 0), (2, 3), (3, 2)}
    // save a bit of time by iterating over pairs
    // This still works if we have equal edge weights because of the way that edge comparison is defined
    for (size_t i = 0; i < edges.size(); i += 2) {
        if (!makesCycle(edges[i].from(), edges[i].to(), mst)) {
            mst[edges[i].from()].push_back(edges[i]);
            mst[edges[i].to()].push_back(edges[i + 1]);
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

    // Start the search at one of the new nodes as we know that if a new cycle is introduced, it will be connected to
    // this node
    return containsCycle(adjList, node1);
}

bool containsCycle(const std::vector<std::list<Edge>>& adjacencyList, const size_t start)
{
    // Could do a short circuit-test here; if there are more edges than nodes, then there is definitely a cycle

    size_t                                numNodes = adjacencyList.size();
    std::vector<bool>                     visited(numNodes, false);
    std::stack<std::pair<size_t, size_t>> stack; // Stack tracking a vertex + its parent

    // Iterate over all our nodes
    // First node does not have a parent so assign something ridiculously big
    visited[start] = true;
    stack.push(std::make_pair(start, -1));
    while (!stack.empty()) {
        // find the node at the top of the stack + remove it
        size_t lastNode = stack.top().first;
        size_t parent   = stack.top().second;
        stack.pop();

        // Iterate over all edges attached to the node
        for (const Edge& edge : adjacencyList[lastNode]) {

            // If this edge attaches to a non visited node, add it to the stack
            if (!visited[edge.to()]) {
                stack.push(std::make_pair(edge.to(), lastNode));
                visited[edge.to()] = true;
            }

            // Otherwise if the edge attaches to a node that is not the parent, we have a cycle
            else if (edge.to() != parent) {
                return true;
            }
        }
    }
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
    size_t order{treeAdjacencyList.size()};
    if (containsCycle(treeAdjacencyList, root)) {
        throw CyclicGraph();
    }

    if (!isSpanning(treeAdjacencyList)) {
        throw NotSpanning();
    }

    // Set our vector of discovered nodes to false, except the one we start at
    std::vector<bool> discovered(order, false);
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
