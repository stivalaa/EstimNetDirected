/*****************************************************************************
 * 
 * File:    undirectedStats.c
 * Author:  Alex Stivala
 * Created: June 2024
 *
 *
 * statistics functions (summing change statistics is verified against these)
 *
 ****************************************************************************/

#include <assert.h>
#include <math.h>
#include "changeStatisticsUndirected.h"
#include "undirectedStats.h"

/*
 * Statistic for FourCycles, number of four-cycles in an undirected graph.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      number of four-cycles in undirected graph g
 */
ulonglong_t FourCycles(const graph_t *g)
{
  uint_t i,l;
  ulonglong_t value = 0;

  assert(!g->is_bipartite);
  assert(!g->is_directed);

  for (i = 1; i < g->num_nodes; i++) {
    for (l = 0; l < i; l++) {
      value += n_choose_2(GET_2PATH_ENTRY(g, i, l));
    }
  }
  assert(value % 2 == 0);
  return value / 2;
}




/*
 * Count number of four-cycles that a particular node u is involved in.
 *
 * Note can also be used as for bipartite networks.
 */
uint_t num_four_cycles_node_SLOW(const graph_t *g, uint_t u)
{
  uint_t v;
  uint_t count = 0;

  /* slow version that iterates over all nodes */
  
  for (v = 0; v < g->num_nodes; v++){
    if (v != u) {
      count += n_choose_2(GET_2PATH_ENTRY(g, u, v));
    }
  }
  return count;
}





/*
 * Statistic for number of 4-cycles at each node raised to a
 * power. The lambda parameter (> 1.0) (mis)used to specify the value
 * 1/lambda as the epxonent. Note this is not the same meaning of
 * lambda as its original use in the "alternating" parameters.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      
 */
double PowerFourCycles(const graph_t *g, double lambda)
{
  uint_t  i;
  double  alpha = 1/lambda;
  ulong_t fourcycle_count = 0;
  double  value = 0;

  assert(!g->is_directed);

  for (i = 0; i < g->num_nodes; i++) {
    fourcycle_count = num_four_cycles_node(g, i);
    uint_t fourcycle_count_SLOW = num_four_cycles_node_SLOW(g, i);
    assert(fourcycle_count_SLOW == fourcycle_count);
    value += pow(fourcycle_count, alpha);
  }
  return value;
}

