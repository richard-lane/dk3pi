#ifndef GRAPHTHEORY_H
#define GRAPHTHEORY_H

#include <list>
#include <vector>

struct CyclicGraph : public std::exception {
    const char* what() const throw() { return "Graph contains a cycle"; }
};

struct NotSpanning: public std::exception {
    const char* what() const throw() { return "Graph not spanning"; }
};

/*
 * Represent an edge by its weight and the nodes that it connects
 *
 * Zero-indexed node numbering
 */
class Edge
{
  public:
    explicit Edge(const size_t fromNode, const size_t toNode, const double weight)
        : _weight(weight), _fromNode(fromNode), _toNode(toNode){};

    /*
     * Comparisons only care about weight
     */
    inline bool operator<(const Edge& other) const { return _weight < other._weight; };
    inline bool operator>(const Edge& other) const { return !(*this < other); };

    /*
     * Equality cares about weight + nodes
     */
    inline bool operator==(const Edge& other) const
    {
        return _weight == other._weight && _fromNode == other.from() && _toNode == other.to();
    };
    inline bool operator!=(const Edge& other) const { return !(*this == other); };

    double getWeight() const { return _weight; };
    double from() const { return _fromNode; };
    double to() const { return _toNode; };

  private:
    double _weight{0.0};
    size_t _fromNode{0};
    size_t _toNode{0};
};

/*
 * Represent an undirected graph by an adjacency list that tracks which nodes are connected to which others
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
    std::vector<std::list<Edge>> getMaxSpanningTree() const;

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
};

/*
 * Check whether adding an edge between two nodes would create a cycle in our graph
 */
bool makesCycle(const size_t node1, const size_t node2, const std::vector<std::list<Edge>>& adjacencyList);

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
 *
 * Only tested on an adirectional graph so it might not work for directed graphs
 *
 * Parent is arbitrary if starting a new search
 */
bool containsCycle(const std::vector<std::list<Edge>>& adjacencyList,
                   const size_t                        node,
                   std::vector<bool>&                  discovered,
                   const size_t                        parent = -1);

/*
 * Check whether the adjacency list representation of a graph is spanning
 *
 * I was tired when i wrote this so it's slow and likely buggy
 */
bool isSpanning(const std::vector<std::list<Edge>>& adjacencyList);

/*
 * Take a tree and find it's in-tree representation
 */
std::vector<std::list<Edge>> inTree(const size_t root, const std::vector<std::list<Edge>>& treeAdjacencyList);

#endif // GRAPHTHEORY_H
