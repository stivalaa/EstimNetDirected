/*****************************************************************************
 * 
 * File:    testStatsSumChangeStatsBipartiteMain.c
 * Author:  Alex Stivala
 * Created: January 2024
 *
 * Test change statistics implementations by comparing sum of change
 * stats for all edges in network to statistic value computed according
 * to the definition of the statistic (implemented in this test module).
 *
 *
 * Usage:  testStatsSumChangeStatsBipartite  <in_edgelistfile> <lambda>
 *
 * Reads graph from Pajek format <in_edgelistfile> and compute stats with
 * weighting parameter <lambda> where real value lambda > 1.
 *
 * Outputs observed statistics value for the statistics, which are computed
 * by summing the change stats over all edges in the data, and verifies
 * that these indeed sum to the statistic value computed directly (in this
 * code).
 *
 * Example:
 * ./testStatsSumChangeStatsBipartite ../../examples/bipartite/simulated/bpnet_A6000_B750_sparse_sim100000000.net 2.0
 *
 *
 ****************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "graph.h"
#include "changeStatisticsBipartiteUndirected.h"
#include "loadGraph.h"

/* Approximate double floating point equality */
#define DOUBLE_APPROX_EQ_TEST(a, b) ( fabs((a) - (b)) <= 1e-08 )



/*****************************************************************************
 *
 * statistics functions (summing change statistics is verified against these)
 *
 ****************************************************************************/


/*
 * Statistic for BipartiteAltKCyclesA, alternating k-cycles for type A
 * nodes (K-Ca in bPNet, XACA in MPNet) defined by eqn (6.12) in:
 *
 *   Wang, P., Sharpe, K., Robins, G. L., & Pattison,
 *   P. E. (2009). Exponential random graph (p∗) models for affiliation
 *   networks. Social Networks, 31(1), 12-25.
 *
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *     lambda - decay value > 1.0
 *
 * Return value:
 *      statistic for g with decay value lambda

 */
static double BipartiteAltKCyclesA(const graph_t *g,double lambda)
{
  uint_t i,j;
  double value = 0;

  assert (lambda > 1.0);

  assert(g->is_bipartite);
  assert(!g->is_directed);
  /* TODO implement */
  return value;
}


/*****************************************************************************
 *
 * main
 *
 ****************************************************************************/

int main(int argc, char *argv[])
{
  uint_t i;
  char *edgelist_filename = NULL;
  FILE *file           = NULL;
  uint_t     num_nodes = 0;
  uint_t     num_A_nodes = 0;
  /*  uint_t     num_P_nodes = 0; */
  graph_t *g         = NULL;
  double stat_value;
  double lambda;
  char  *endptr; /* for strtod() */
 
  srand(time(NULL));

  if (argc != 3) {
    fprintf(stderr, "Usage: %s <edgelist_file> <lambda>\n", argv[0]);
    exit(1);
  }
  edgelist_filename = argv[1];
  lambda = strtod(argv[2], &endptr);
  if (lambda <= 1.0) {
    fprintf(stderr, "lambda value %g is not > 1.0\n", lambda);
    return -1;
  }

  if (!(file = fopen(edgelist_filename, "r"))) {
    fprintf(stderr, "error opening file %s (%s)\n",
            edgelist_filename, strerror(errno));
    return -1;
  }
  
  get_num_vertices_from_bipartite_pajek_file(file,
					     &num_nodes,
					     &num_A_nodes);/* closes file */
  
  g = allocate_graph(num_nodes, FALSE/*undirected*/, TRUE/*bipartite*/,
		     num_A_nodes);

  if (!(file = fopen(edgelist_filename, "r"))) {
    fprintf(stderr, "error opening file %s (%s)\n",
            edgelist_filename, strerror(errno));
    return -1;
  }



#define NUM_FUNCS 1
  uint_t n_total = NUM_FUNCS;
  static double lambda_values[NUM_FUNCS];
  double obs_stats[NUM_FUNCS];
  static double theta[NUM_FUNCS]; /* init to zero, unused */
  change_stats_func_t *change_stats_funcs[NUM_FUNCS];


  change_stats_funcs[0] = &changeBipartiteAltKCyclesA;
  lambda_values[0]      = lambda;

  for (i = 0; i < NUM_FUNCS; i++) {
    obs_stats[i] = 0;
  }
  g = load_graph_from_arclist_file(file, g, TRUE,
                                   n_total, 0, 0, 0, change_stats_funcs,
                                   lambda_values, NULL,
                                   NULL, NULL, NULL, NULL,
                                   NULL, obs_stats, theta);
  for (i = 0; i < NUM_FUNCS; i++) {
    printf("%g ", obs_stats[i]);
  }
  printf("\n");

  stat_value= BipartiteAltKCyclesA(g, lambda_values[0]);
  assert(DOUBLE_APPROX_EQ(stat_value,  obs_stats[0]));


  free_graph(g);

  
  exit(0);
}
