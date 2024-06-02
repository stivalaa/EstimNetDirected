#ifndef BIPARTITESTATS_H
#define BIPARTITESTATS_H
/*****************************************************************************
 * 
 * File:    bipartiteStats.h
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

double FourCyclesA(const graph_t *g);
double FourCyclesB(const graph_t *g);
double BipartiteAltKCyclesA(const graph_t *g, double lambda);
double BipartiteAltKCyclesB(const graph_t *g, double lambda);
double BipartiteAltKCyclesA_SLOW(const graph_t *g, double lambda);
double BipartiteAltKCyclesB_SLOW(const graph_t *g, double lambda);
double BipartiteAltK4CyclesA_SLOW(const graph_t *g, double lambda);
double PowerFourCyclesA(const graph_t *g, double lambda);
double PowerFourCyclesB(const graph_t *g, double lambda);
double PowerFourCycles(const graph_t *g, double lambda);


#endif /* BIPARTITESTATS_H */
