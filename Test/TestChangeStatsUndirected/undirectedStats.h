#ifndef UNDIRECTEDSTATS_H
#define UNDIRECTEDSTATS_H
/*****************************************************************************
 * 
 * File:    undirectedStats.h
 * Author:  Alex Stivala
 * Created: June 2024
 *
 *
 * statistics functions (summing change statistics is verified against these)
 *
 *
 ****************************************************************************/
#include "utils.h"
#include "graph.h"


/* Approximate double floating point equality */
#define DOUBLE_APPROX_EQ_TEST(a, b) ( fabs((a) - (b)) <= 1e-06 )


ulonglong_t FourCycles(const graph_t *g);
uint_t num_four_cycles_node_SLOW(const graph_t *g, uint_t u);
double PowerFourCycles(const graph_t *g, double lambda);

#endif /* UNDIRECTEDSTATS_H */
