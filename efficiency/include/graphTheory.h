#ifndef GRAPHTHEORY_H
#define GRAPHTHEORY_H

#include <list>
#include <vector>

/*
 * Represent an edge by its weight and the node that it connects to
 *
 * Zero-indexed node numbering
 *
 * This class might have way more information in it than we need
 */
class Edge
{
  public:
    explicit Edge(const size_t node, const double weight);

    bool operator==(const Edge& other) const;
    bool operator!=(const Edge& other) const;
    bool operator<(const Edge& other) const;
    bool operator>(const Edge& other) const;

    double getWeight() const { return _weight; };
    double node() const { return _node; };

  private:
    const double _weight{0.0};
    const size_t _node{0};
};

/*
 * Represent a graph by an adjacency list that tracks which nodes are connected to which others
 */
class Graph
{
  public:
    explicit Graph(const size_t numNodes);

    /*
     * Add an edge connecting the two nodes indexed by node1 and node2
     *
     * Maybe this should take an Edge object but nah
     */
    void addEdge(const size_t node1, const size_t node2, const double weight);

    /*
     * Return the adjacency list representation of the minimum spanning tree of this graph
     */
    std::vector<std::list<Edge>> getMST() const;

    std::vector<std::list<Edge>> getAdjacencyList() const { return _adjacencyList; };

  private:
    size_t _numNodes{0};

    /*
     * Represent our graph by an array of linked lists
     *
     * if there is only edge, connecting nodes 1 and 3 (indexed from 0), then this will look like:
     *     {{}, {Edge(3, w13)}, {}, {Edge(1, w13), {}, {}, {}...}}
     */
    std::vector<std::list<Edge>> _adjacencyList{};

    /*
     * Check whether adding an edge between two nodes would create a cycle in our graph
     *
     */
    bool _makesCycle(const size_t node1, const size_t node2) const;
};

/*
 * Search a graph depth-first for cycles
 *
 * Must provide a representation of the graph as an adjacency list, the node to start the search from, a std::vector
 * representing whether this node has been visited before and the index of the node's parent.
 *
 * To start a new search, discovered should be a vector of "false"
 *
 * Starts at the provided node and searches each node connected to this one in order recursively
 *
 * Will not detect cycles that are not connected to the node initially passed in
 *
 * Also is a recursive function so could cause a stack overflow if the graph is huge
 */
bool containsCycle(const std::vector<std::list<Edge>>& adjacencyList,
                   const size_t                        node,
                   std::vector<bool>&                  discovered,
                   const size_t                        parent);

/*
 * Check whether the adjacency list representation of a graph is spanning
 *
 * I was tired when i wrote this so it's slow and likely buggy
 */
bool isSpanning(const std::vector<std::list<Edge>>& adjacencyList);

#endif // GRAPHTHEORY_H
