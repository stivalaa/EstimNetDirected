/*****************************************************************************
 * 
 * File:    testStatsSumChangeStatsAttrBipartiteMain.c
 * Author:  Alex Stivala
 * Created: January 2025
 *
 * Test the bipartite attribute change statistics by comparing
 * directly calculated statistics with values obtained by summing
 * change statistics over all edges in network.
 *
 * Usage:  testStatsSumChangeStatsAttrBipartite  <in_edgelistfile> <binattr_file> 
 *
 * Reads graph from Pajek format <in_edgelistfile> and binary
 * attributes from <binattr_file>.
 *
 * Outputs observed statistics value for the statistics, which are computed
 * by summing the change stats over all edges in the data, and verifies
 * that these indeed sum to the statistic value computed directly (in this
 * code). Outputs the observed statistics for external validation
 * (by scripts with hardcoded manually checked values and/or comparison
 * to results from statnet for example)
 *
 * Example:
 *   testStatsSumChangeStatsAttrBipartite   ../../examples/bipartite/simulated/bpnet_A12000_B4000_attrs_sim830000000.ne  ../../examples/bipartite/simulation/binattr_all.txt
 *
 ****************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "changeStatisticsBipartiteUndirected.h"
#include "loadGraph.h"
#include "attrbipartiteStats.h"


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
  char *binattr_filename;
  double stat_value;
 
  srand(time(NULL));

  if (argc != 3) {
    fprintf(stderr, "Usage: %s <inedgelist_file> <binattr_file>\n", argv[0]);
    exit(1);
  }
  edgelist_filename = argv[1];
  binattr_filename = argv[2];


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

  if (load_attributes(g, binattr_filename, NULL, NULL, NULL) != 0) {
    fprintf(stderr, "ERRROR: load node attributes failed\n");
    exit(1);
  }

  if (!(file = fopen(edgelist_filename, "r"))) {
    fprintf(stderr, "error opening file %s (%s)\n",
            edgelist_filename, strerror(errno));
    return -1;
  }

  
  /* hardcoding indices of attributes to match input files*/
  /* binattr_all.txt: binattrA binattrP binattrAP */
  const uint_t binattrA_index = 0;
  const uint_t binattrP_index = 1;
  const uint_t binattrAP_index = 2;



#define NUM_FUNCS 3
  uint_t n_total = NUM_FUNCS, n_attr = NUM_FUNCS;
  uint_t attr_indices[NUM_FUNCS];
  static double lambda_values[NUM_FUNCS]; /* init to zero, unused */
  static double exponent_values[NUM_FUNCS]; /* init to zero, unused */
  double obs_stats[NUM_FUNCS];
  static double theta[NUM_FUNCS]; /* init to zero, unused */
  attr_change_stats_func_t *attr_change_stats_funcs[NUM_FUNCS];


  
  attr_change_stats_funcs[0] = &changeBipartiteExactlyOneNeighbourA;
  attr_indices[0]            = binattrP_index;

  attr_change_stats_funcs[1] = &changeBipartiteExactlyOneNeighbourB;
  attr_indices[1]            = binattrA_index;

  attr_change_stats_funcs[2] = &changeBipartiteTwoPathExactlyOneNeighbourA;
  attr_indices[2]            = binattrP_index;

  
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

  stat_value= BipartiteExactlyOneNeighbourA(g, attr_indices[0]);
  assert(DOUBLE_APPROX_EQ(stat_value,  obs_stats[0]));

  stat_value= BipartiteExactlyOneNeighbourB(g, attr_indices[1]);
  assert(DOUBLE_APPROX_EQ(stat_value,  obs_stats[1]));

  stat_value= BipartiteTwoPathExactlyOneNeighbourA(g, attr_indices[2]);
  assert(DOUBLE_APPROX_EQ(stat_value,  obs_stats[2]));

  free_graph(g);

  exit(0);
}
