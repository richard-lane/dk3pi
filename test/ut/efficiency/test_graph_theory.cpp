#ifndef TEST_GRAPH_THEORY
#define TEST_GRAPH_THEORY

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "graphTheory.h"

/*
 * Test Edge class equality
 */
BOOST_AUTO_TEST_CASE(test_edge_equality)
{
    // Check that edges with the same weight and nodes are equal
    BOOST_CHECK(Edge(0, 1, 4.0) == Edge(0, 1, 4.0));

    // Check edges with the same weight but different nodes are not equal
    BOOST_CHECK(Edge(0, 1, 4.0) != Edge(1, 0, 4.0));
    BOOST_CHECK(Edge(0, 1, 4.0) != Edge(2, 3, 4.0));

    // Check edges with different weight but the same nodes are not equal
    BOOST_CHECK(Edge(0, 1, 4.0) != Edge(1, 0, 8.0));
}

/*
 * Test Edge class comparison operators
 */
BOOST_AUTO_TEST_CASE(test_edge_comparison)
{
    BOOST_CHECK(Edge(0, 1, 4.0) < Edge(0, 1, 5.0));
    BOOST_CHECK(Edge(2, 3, 6.0) > Edge(0, 1, 5.0));
}

/*
 * Test that we can correctly identify whether adding an edge to a graph will make a cycle
 */
BOOST_AUTO_TEST_CASE(test_makes_cycle)
{
    // Create a graph thing that doesn't contain a cycle
    const double                       wt{1.0};
    const std::vector<std::list<Edge>> graph = {
        {Edge(0, 1, wt), Edge(0, 2, wt)}, {Edge(1, 0, wt), Edge(1, 4, wt)}, {Edge(2, 0, wt)}, {}, {Edge(4, 1, wt)}};

    BOOST_CHECK(makesCycle(2, 4, graph));
    BOOST_CHECK(!makesCycle(3, 4, graph));
}

/*
 * Test we correctly detect a cycle in a graph
 */
BOOST_AUTO_TEST_CASE(test_find_cycles)
{
    // Create a graph thing that doesn't contain a cycle
    const double                       wt{1.0};
    const std::vector<std::list<Edge>> noCycle = {{Edge(0, 1, wt), Edge(0, 2, wt)},
                                                  {Edge(1, 0, wt), Edge(1, 3, wt), Edge(1, 4, wt)},
                                                  {Edge(2, 0, wt)},
                                                  {Edge(3, 1, wt)},
                                                  {Edge(4, 1, wt)}};

    const std::vector<std::list<Edge>> yesCycle = {{Edge(0, 1, wt), Edge(0, 2, wt)},
                                                   {Edge(1, 0, wt), Edge(1, 3, wt), Edge(1, 4, wt)},
                                                   {Edge(2, 0, wt), Edge(2, 4, wt)},
                                                   {Edge(3, 1, wt)},
                                                   {Edge(4, 1, wt), Edge(4, 2, wt)}};

    std::vector<bool> discovered(5, false);
    BOOST_CHECK(!containsCycle(noCycle, 0, discovered, -1));

    discovered = std::vector<bool>(5, false);
    BOOST_CHECK(containsCycle(yesCycle, 0, discovered, -1));
}

/*
 * Test that we correctly detemine whether a graph is spanning
 */
BOOST_AUTO_TEST_CASE(test_spanning)
{
    const double                       wt{1.0};
    const std::vector<std::list<Edge>> connected = {{Edge(0, 1, wt), Edge(0, 3, wt)},
                                                    {Edge(1, 0, wt)},
                                                    {Edge(2, 3, wt), Edge(2, 4, wt)},
                                                    {Edge(3, 0, wt), Edge(3, 2, wt)},
                                                    {Edge(4, 2, wt)}};

    const std::vector<std::list<Edge>> disconnected = {
        {Edge(0, 0, wt), Edge(0, 1, wt)}, {Edge(1, 2, wt)}, {Edge(2, 1, wt)}, {Edge(3, 4, wt)}, {Edge(4, 3, wt)}};

    BOOST_CHECK(isSpanning(connected));
    BOOST_CHECK(!isSpanning(disconnected));
}

/*
 * Test that we correctly identify the maximum spanning tree
 */
BOOST_AUTO_TEST_CASE(test_kruskal)
{
    // Our edges will get added in decreasing order of weight
    const std::vector<std::list<Edge>> expectedMaxSpanningTree = {{Edge(0, 1, 8)},
                                                                  {Edge(1, 3, 11), Edge(1, 0, 8)},
                                                                  {Edge(2, 3, 15), Edge(2, 4, 10)},
                                                                  {Edge(3, 2, 15), Edge(3, 1, 11)},
                                                                  {Edge(4, 2, 10)}};
    // Create and draw edges on our graph
    Graph graph(5);
    graph.addEdge(0, 1, 8);
    graph.addEdge(0, 2, 5);
    graph.addEdge(1, 2, 9);
    graph.addEdge(1, 3, 11);
    graph.addEdge(2, 3, 15);
    graph.addEdge(2, 4, 10);
    graph.addEdge(3, 4, 7);

    BOOST_CHECK(graph.getMaxSpanningTree() == expectedMaxSpanningTree);
}

#endif // TEST_GRAPH_THEORY
