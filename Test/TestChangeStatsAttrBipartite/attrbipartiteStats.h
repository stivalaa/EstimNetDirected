#ifndef ATTRBIPARTITESTATS_H
#define ATTRBIPARTITESTATS_H

/*****************************************************************************
 * 
 * File:    attrbipartiteStats.h
 * Author:  Alex Stivala
 * Created: January 2025
 *
 *
 * statistics functions (summing change statistics is verified against these)
 *
 ****************************************************************************/

#include "utils.h"
#include "graph.h"


/* Approximate double floating point equality */
#define DOUBLE_APPROX_EQ_TEST(a, b) ( fabs((a) - (b)) <= 1e-08 )

double BipartiteExactlyOneNeighbourA(const graph_t *g, uint_t a);
double BipartiteExactlyOneNeighbourB(const graph_t *g, uint_t a);
double BipartiteTwoPathExactlyOneNeighbourA(const graph_t *g, uint_t a);
double BipartiteTwoPathExactlyOneNeighbourB(const graph_t *g, uint_t a);

#endif /* ATTRBIPARTITESTATS_H */
