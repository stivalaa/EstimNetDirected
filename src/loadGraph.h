#ifndef LOADGRAPH_H
#define LOADGRAPH_H
/*****************************************************************************
 * File:    loadGraph.h
 * Author:  Alex Stivala
 * Created: October 2017
 *
 * Load digraph or graph from Pajek format arc list file and optionally compute
 * statistics corresponding to ERGM parameters.
 *
 ****************************************************************************/

#include "graph.h"
#include "changeStatisticsGeneral.h"

graph_t *load_graph_from_arclist_file(FILE *pajek_file, graph_t *g,
                                          bool computeStats,
                                          uint_t n, uint_t n_attr,
                                          uint_t n_dyadic,
                                          uint_t n_attr_interaction,
                                          change_stats_func_t
                                                     *change_stats_funcs[],
                                          double lambda_values[],
                                          attr_change_stats_func_t
                                          *attr_change_stats_funcs[],
                                          dyadic_change_stats_func_t
                                          *dyadic_change_stats_funcs[],
                                          attr_interaction_change_stats_func_t
                                         *attr_interaction_change_stats_funcs[],
                                          uint_t attr_indices[],
                                          double exponent_values[],
                                          uint_pair_t
                                          attr_interaction_pair_indices[],
                                          double addChangeStats[],
                                          double theta[]);


#endif /* LOADGRAPH_H */

