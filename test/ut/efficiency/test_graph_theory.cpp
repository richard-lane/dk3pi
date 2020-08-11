#ifndef TEST_GRAPH_THEORY
#define TEST_GRAPH_THEORY

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "graphTheory.h"

/*
 * Test we correctly detect a cycle in a graph
 */
BOOST_AUTO_TEST_CASE(test_find_cycles)
{
    // Create a graph thing that doesn't contain a cycle
    const double                       wt{1.0};
    const std::vector<std::list<Edge>> noCycle = {{Edge(1, wt), Edge(2, wt)},
                                                  {Edge(0, wt), Edge(3, wt), Edge(4, wt)},
                                                  {Edge(0, wt)},
                                                  {Edge(1, wt)},
                                                  {Edge(1, wt)}};

    const std::vector<std::list<Edge>> yesCycle = {{Edge(1, wt), Edge(2, wt)},
                                                   {Edge(0, wt), Edge(3, wt), Edge(4, wt)},
                                                   {Edge(0, wt), Edge(4, wt)},
                                                   {Edge(1, wt)},
                                                   {Edge(1, wt), Edge(2, wt)}};

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
    const std::vector<std::list<Edge>> connected = {
        {Edge(1, wt), Edge(3, wt)}, {Edge(0, wt)}, {Edge(3, wt), Edge(4, wt)}, {Edge(0, wt), Edge(2, wt)}, {Edge(2, wt)}};

    const std::vector<std::list<Edge>> disconnected = {
        {Edge(0, wt), Edge(1, wt)}, {Edge(2, wt)}, {Edge(1, wt)}, {Edge(4, wt)}, {Edge(3, wt)}};

    BOOST_CHECK(isSpanning(connected));
    BOOST_CHECK(!isSpanning(disconnected));
}

#endif // TEST_GRAPH_THEORY
