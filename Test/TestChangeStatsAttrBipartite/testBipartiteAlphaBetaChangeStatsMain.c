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
 * Outputs observed statistics value for the statistics, which are computed
 * by summing the change stats over all edges in the data, and verifies
 * that these indeed sum to the statistic value computed directly (in this
 * code). Outputs the observed statistics for external validation
 * (by scripts with hardcoded manually checked values and/or comparison
 * to results from statnet for example)
 *
 * Example:
 * ./testBipartiteAlphaBetaChangeStats ../../examples/bipartite/simulated/bpnet_A12000_B4000_attrs_sim830000000.net ../../examples/bipartite/simulation/catattr_all.txt
 *
 * b1nodematch and b2nodematch (statnet ergm names) are defined in:
 *
 *    Bomiriya, R. P. (2014). Topics in exponential random graph
 *    modeling. (Doctoral dissertation, Pennsylvania State University).
 *    https://etda.libraries.psu.edu/files/final_submissions/9857
 *
 *    Bomiriya, R. P., Kuvelkar, A. R., Hunter, D. R., & Triebel,
 *    S. (2023). Modeling Homophily in Exponential-Family Random Graph
 *    Models for Bipartite Networks. arXiv preprint
 *    arXiv:2312.05673. https://arxiv.org/abs/2312.05673
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

/*****************************************************************************
 *
 * statistics functions (summing change statistics is verified against these)
 *
 ****************************************************************************/

/*
 * Statistic for Bipartite node-centered (alpha-based) homophily
 * for type A node (b1nodematch(alpha) statnet ergm term)
 *
 * alpha is the exponent in the range [0, 1]
 *
 * b1nodematch and b2nodematch (statnet ergm names) are defined in:
 *
 *  Bomiriya, R. P. (2014). Topics in exponential random graph
 *  modeling. (Doctoral dissertation, Pennsylvania State University).
 *  https://etda.libraries.psu.edu/files/final_submissions/9857
 *
 *  Bomiriya, R. P., Kuvelkar, A. R., Hunter, D. R., & Triebel,
 *  S. (2023). Modeling Homophily in Exponential-Family Random Graph
 *  Models for Bipartite Networks. arXiv preprint
 *  arXiv:2312.05673. https://arxiv.org/abs/2312.05673
 *
 * This statistic is defined by equation (6) in Bomiriya et al. (2023)
 *
 * Parameters:
 *     g     - undirected bipartite graph
 *     a     - index of categorical attribute
 *     alpha - exponent in range [0, 1]
 *
 * Return value:
 *     node-centered (alpha-based) homophily statistic for g with cat attr a

 */
static double BipartiteNodematchAlphaA(const graph_t *g, uint_t a,
				       double alpha)
{
  uint_t i,j;
  double value = 0;

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = 0; i < g->num_A_nodes; i++) {
    for (j = 0; j < i; j++) { /* do not double-count (i,j) two-paths */
      if (j != i &&
	  g->catattr[a][i] != CAT_NA &&
	  g->catattr[a][j] != CAT_NA &&
	  g->catattr[a][i] == g->catattr[a][j]) {
	/* Note pow0 defines pow0(0, 0) = 0
           as per Bomiryia et al. (2023) [see p. 7 after eqn (7)] */
	value += pow0(GET_A2PATH_ENTRY(g, i, j), alpha);
      }
    }
  }
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
  struct timeval start_timeval, end_timeval, elapsed_timeval;
  int    etime;
  char *catattr_filename;
  double stat_value;
 
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
  exponent_values[0]         = 1;
  
  attr_change_stats_funcs[1] = &changeBipartiteNodematchBetaA;
  attr_indices[1]            = catattrA_index;
  exponent_values[1]         = 1;

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

  stat_value= BipartiteNodematchAlphaA(g, attr_indices[0], exponent_values[0]);
  /*printf("%g\n", stat_value);*/
  assert(DOUBLE_APPROX_EQ(stat_value,  obs_stats[0]));
  
  free_graph(g);
  exit(0);
}
