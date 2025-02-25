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
 * Statistic for FourCycles, number of four-cycles in an undirected graph.
 * This version computes it by summing number of four-cycles at each node
 * compuited with num_four_cycles_node().
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      number of four-cycles in undirected graph g
 */
ulonglong_t FourCycles_sum_by_node(const graph_t *g)
{
  uint_t i;
  ulonglong_t value = 0;

  assert(!g->is_bipartite);
  assert(!g->is_directed);

  for (i = 0; i < g->num_nodes; i++) {
    value += num_four_cycles_node(g, i);
  }
  /* Each four-cycle is counted 4 times, once for each node in it */
  assert(value % 4 == 0);
  return value / 4;
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
 * For the definition of the statistic, see:
 *
 *   Stivala, A., Wang, P., & Lomi, A. (2025). Improving
 *   exponential-family random graph models for bipartite networks.
 *   arXiv preprint arXiv:2502.01892. https://arxiv.org/abs/2502.01892
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      FourCyclesNodePower statistic
 *
 */
double PowerFourCycles(const graph_t *g, double lambda)
{
  uint_t  i;
  double  alpha = 1/lambda;
  ulong_t fourcycle_count = 0;
  ulonglong_t fourcycle_count_sum = 0;
  ulonglong_t num_fourcycles;
  double  value = 0;

  assert(!g->is_directed);

  for (i = 0; i < g->num_nodes; i++) {
    fourcycle_count = num_four_cycles_node(g, i);
    uint_t fourcycle_count_SLOW = num_four_cycles_node_SLOW(g, i);
    assert(fourcycle_count_SLOW == fourcycle_count);
    fourcycle_count_sum += fourcycle_count;
    value += pow(fourcycle_count, alpha);
  }
  /* Each four-cycle countains 4 nodes so is counted 4 times */
  assert(fourcycle_count_sum % 4 == 0);
  //fprintf(stderr, "fourcycle_count_sum = %llu\n", fourcycle_count_sum);
  num_fourcycles = FourCycles(g);
  //fprintf(stderr, "num_fourcycles =      %llu\n", num_fourcycles);
  assert(fourcycle_count_sum / 4 == num_fourcycles);
  return value;
}

