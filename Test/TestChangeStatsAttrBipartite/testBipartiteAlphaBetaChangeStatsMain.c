/*****************************************************************************
 * 
 * File:    testBipartiteAlphaBetaChangeStatsMain.c
 * Author:  Alex Stivala
 * Created: December 2023
 *
 * Test the BipartiteNodeMatchAlpha[AB] and BipartiteNodeMatchBeta[AB]
 * (statnet ergm b1nodematch and b2nodematch) change statistics.
 *
 *
 * Usage:  testBipartiteAlphaBetaChangeStats  <in_edgelistfile>  <catattr_file>
 *
 * Reads graph from Pajek format <in_edgelistfile> and attributes from
 *  <catattr_file>
 *
 * Outputs observed statistics value for the statistics, which are compute
 * by summing the change stats over all edges in the data, and verifies
 * that these indeed sum to the statistic value computed directly (in this
 * code). Outputs the observed statistics for external validation
 * (by scripts with hardcoded manually checked values and/or comparison
 * to results from statnet for example)
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


int main(int argc, char *argv[]) 
{
  uint_t i,j;
  char *edgelist_filename = NULL;
  FILE *file           = NULL;
  uint_t     num_nodes = 0;
  uint_t     num_A_nodes = 0;
  /*  uint_t     num_P_nodes = 0; */
  graph_t *g         = NULL;
  struct timeval start_timeval, end_timeval, elapsed_timeval;
  int    etime;
  char *catattr_filename;
 
  srand(time(NULL));

  if (argc != 3) {
    fprintf(stderr, "Usage: %s <inedgelist_file> <catattr_file>\n", argv[0]);
    exit(1);
  }
  edgelist_filename = argv[1];
  catattr_filename = argv[2];

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

  if (load_attributes(g, NULL, catattr_filename, NULL,
		      NULL) != 0) {
    fprintf(stderr, "ERRROR: load node attributes failed\n");
    exit(1);
  }

  if (!(file = fopen(edgelist_filename, "r"))) {
    fprintf(stderr, "error opening file %s (%s)\n", 
            edgelist_filename, strerror(errno));
    return -1;
  }



  
  /* hardcoding indices of attributes to match input files*/
  /* catattr_all.txt: catattrA catattrP catattrAP */
  const uint_t catattrA_index = 0;
  const uint_t catattrP_index = 1;
  const uint_t catattrAP_index = 2;

#define NUM_FUNCS 2
  uint_t n_total = NUM_FUNCS, n_attr = NUM_FUNCS;
  uint_t attr_indices[NUM_FUNCS];
  static double lambda_values[NUM_FUNCS]; /* init to zero, unused */
  double exponent_values[NUM_FUNCS], obs_stats[NUM_FUNCS];
  static double theta[NUM_FUNCS]; /* init to zero, unused */
  attr_change_stats_func_t *attr_change_stats_funcs[NUM_FUNCS];

  
  attr_change_stats_funcs[0] = &changeBipartiteNodematchAlphaA;
  attr_indices[0]            = catattrA_index;
  exponent_values[0]         = 0.1;
  
  attr_change_stats_funcs[1] = &changeBipartiteNodematchBetaA;
  attr_indices[1]            = catattrA_index;
  exponent_values[1]         = 0.1;

  for (i = 0; i < NUM_FUNCS; i++) {
    obs_stats[i] = 0;
  }
  g = load_graph_from_arclist_file(file, g, TRUE,
                                   n_total, n_attr, 0, 0, NULL,
                                   lambda_values, attr_change_stats_funcs,
                                   NULL, NULL, attr_indices, exponent_values,
                                   NULL, obs_stats, theta);
  for (i = 0; i < NUM_FUNCS; i++) {
    printf("%g ", obs_stats[i]);
  }
  printf("\n");
  
  free_graph(g);
  exit(0);
}
