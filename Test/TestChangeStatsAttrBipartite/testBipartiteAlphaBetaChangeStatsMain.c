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
 * Usage:  testBipartiteAlphaBetaChangeStats  <in_edgelistfile> <conattr_file> <catattr_file> <exponent>
 *
 * Reads graph from Pajek format <in_edgelistfile> and continuous
 *  attributes from <conattR_file> and categorical attributes from
 *  <catattr_file> and compute stats with alpha or beta value
 *  <exponent> floating point in [0,1]
 *
 * Outputs observed statistics value for the statistics, which are computed
 * by summing the change stats over all edges in the data, and verifies
 * that these indeed sum to the statistic value computed directly (in this
 * code). Outputs the observed statistics for external validation
 * (by scripts with hardcoded manually checked values and/or comparison
 * to results from statnet for example)
 *
 * Example:
 * ./testBipartiteAlphaBetaChangeStats ../../examples/bipartite/simulated/bpnet_A12000_B4000_attrs_sim830000000.net ../../examples/bipartite/simulation/conattr_all.txt  ../../examples/bipartite/simulation/catattr_all.txt 0.1
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

/* Approximate double floating point equality */
#define DOUBLE_APPROX_EQ_TEST(a, b) ( fabs((a) - (b)) <= 1e-08 )

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

  assert (alpha >= 0 && alpha <= 1);

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = 0; i < g->num_A_nodes; i++) {
    assert(bipartite_node_mode(g, i) == MODE_A);    
    for (j = 0; j < i; j++) { /* do not double-count (i,j) two-paths */
      assert(bipartite_node_mode(g, j) == MODE_A);      
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


/*
 * Statistic for Bipartite node-centered (alpha-based) homophily
 * for type B node (b2nodematch(alpha) statnet ergm term)
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
static double BipartiteNodematchAlphaB(const graph_t *g, uint_t a,
				       double alpha)
{
  uint_t i,j;
  double value = 0;

  assert (alpha >= 0 && alpha <= 1);

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = g->num_A_nodes; i < g->num_nodes; i++) {
    assert(bipartite_node_mode(g, i) == MODE_B);
    for (j = g->num_A_nodes; j < i; j++) { /* do not double-count (i,j) two-paths */
      assert(bipartite_node_mode(g, j) == MODE_B);      
      if (j != i &&
	  g->catattr[a][i] != CAT_NA &&
	  g->catattr[a][j] != CAT_NA &&
	  g->catattr[a][i] == g->catattr[a][j]) {
	/* Note pow0 defines pow0(0, 0) = 0
           as per Bomiryia et al. (2023) [see p. 7 after eqn (7)] */
	value += pow0(GET_B2PATH_ENTRY(g, i, j), alpha);
      }
    }
  }
  return value;
}


/*
 * Statistic for Bipartite edge-centered (beta-based) homophily
 * for type A node (b1nodematch(beta) statnet ergm term)
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
 * This statistic is defined by equation (7) in Bomiriya et al. (2023)
 *
 * Parameters:
 *     g     - undirected bipartite graph
 *     a     - index of categorical attribute
 *     alpha - exponent in range [0, 1]
 *
 * Return value:
 *     edge-centered (beta-based) homophily statistic for g with cat attr a

 */
static double BipartiteNodematchBetaA(const graph_t *g, uint_t a,
				      double beta)
{
  uint_t i,j,k,l,m;
  double value = 0;

  assert (beta >= 0 && beta <= 1);
  
  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = 0; i < g->num_A_nodes; i++) {
    assert(bipartite_node_mode(g, i) == MODE_A);
    for (l = 0; l < g->degree[i]; l++) {
      k = g->edgelist[i][l]; /* k iterates over neighbours of i */
      assert(bipartite_node_mode(g, k) == MODE_B);
      uint_t u = 0; /* number of edges to k from nodes (not i) matching i */
      for (m = 0; m < g->degree[k]; m++) {
	j = g->edgelist[k][m]; /* j iterates over neighbours of k */
	assert(bipartite_node_mode(g, j) == MODE_A);
	if (j != i &&
	    g->catattr[a][i] != CAT_NA &&
	    g->catattr[a][j] != CAT_NA &&
	    g->catattr[a][i] == g->catattr[a][j]) {
	  u++;
	}
      }
      /* Note pow0 defines pow0(0, 0) = 0
	 as per Bomiryia et al. (2023) [see p. 7 after eqn (7)] */
      value += pow0(u, beta);
    }
  }
  value /= 2;
  return value;
}

/*
 * Statistic for Bipartite edge-centered (beta-based) homophily
 * for type B node (b2nodematch(beta) statnet ergm term)
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
 * This statistic is defined by equation (7) in Bomiriya et al. (2023)
 *
 * Parameters:
 *     g     - undirected bipartite graph
 *     a     - index of categorical attribute
 *     alpha - exponent in range [0, 1]
 *
 * Return value:
 *     edge-centered (beta-based) homophily statistic for g with cat attr a

 */
static double BipartiteNodematchBetaB(const graph_t *g, uint_t a,
				      double beta)
{
  uint_t i,j,k,l,m;
  double value = 0;

  assert (beta >= 0 && beta <= 1);
  
  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = g->num_A_nodes; i < g->num_nodes; i++) {
    assert(bipartite_node_mode(g, i) == MODE_B);
    for (l = 0; l < g->degree[i]; l++) {
      k = g->edgelist[i][l]; /* k iterates over neighbours of i */
      assert(bipartite_node_mode(g, k) == MODE_A);
      uint_t u = 0; /* number of edges to k from nodes (not i) matching i */
      for (m = 0; m < g->degree[k]; m++) {
	j = g->edgelist[k][m]; /* j iterates over neighbours of k */
	assert(bipartite_node_mode(g, j) == MODE_B);
	if (j != i &&
	    g->catattr[a][i] != CAT_NA &&
	    g->catattr[a][j] != CAT_NA &&
	    g->catattr[a][i] == g->catattr[a][j]) {
	  u++;
	}
      }
      /* Note pow0 defines pow0(0, 0) = 0
	 as per Bomiryia et al. (2023) [see p. 7 after eqn (7)] */
      value += pow0(u, beta);
    }
  }
  value /= 2;
  return value;
}



/*
 * Statistic for Bipartite edge-centered (beta-based)
 * continuous absolute difference (heterophily on continuous
 * attribute) for type A node.
 *
 * beta is the exponent in the range [0, 1]
 *
 * There is no statnet equivalent for this term, but it is based on
 * b1nodematch and b2nodematch (statnet ergm names) as defined in:
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
 * But instead of counting two-paths between matching nodes, it sums
 * the absolute differences between the continuous attributes of nodes
 * connected by two-paths (see changeBipartiteNodeMatchBeta() for the
 * original implementation for matching categorical attributres from which
 * this is derived).
 *
 * This is experimental and not currently used as we do no have
 * a change statistic derived for this statistic, and also it
 * is not clear this statistic is useful: BipartiteNodeMatchBeta
 * makes sense as raising the count of matching noes to a power
 * in [0, 1] down-weights the contribution of each additional
 * matching node, however it is not clear that doing this with a
 * sum of differences in continuous attributes has any sensible
 * interpretation.
 */
static double BipartiteDiffBetaA(graph_t *g, uint_t a, double beta)
{
  uint_t i,j,k,l,m;
  double value = 0;

  assert (beta >= 0 && beta <= 1);

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = 0; i < g->num_A_nodes; i++) {
    assert(bipartite_node_mode(g, i) == MODE_A);
    for (l = 0; l < g->degree[i]; l++) {
      k = g->edgelist[i][l]; /* k iterates over neighbours of i */
      assert(bipartite_node_mode(g, k) == MODE_B);
      double u = 0; /* accumulate sum of abs diff of continuous attributes on
                       connected nodes i -- j where j != i */
      for (m = 0; m < g->degree[k]; m++) {
	j = g->edgelist[k][m]; /* j iterates over neighbours of k */
	assert(bipartite_node_mode(g, j) == MODE_A);
	if (j != i &&
	    !isnan(g->contattr[a][i]) && !isnan(g->contattr[a][j])) {
	  u += fabs(g->contattr[a][i] - g->contattr[a][j]);
	}
      }
      value += pow(u, beta);
    }
  }
  value /= 2;
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
  char *conattr_filename;
  char *catattr_filename;
  double stat_value;
  double exponent;
  char        *endptr; /* for strtod() */
 
  srand(time(NULL));

  if (argc != 5) {
    fprintf(stderr, "Usage: %s <inedgelist_file> <conattr_file> <catattr_file> <exponent>\n", argv[0]);
    exit(1);
  }
  edgelist_filename = argv[1];
  conattr_filename = argv[2];
  catattr_filename = argv[3];
  exponent = strtod(argv[4], &endptr);
  if (exponent < 0.0 || exponent > 1.0) {
    fprintf(stderr, "exponent %g is not in [0, 1]\n", exponent);
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

  if (load_attributes(g, NULL, catattr_filename, conattr_filename,
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

  /* hardcoding indices of attributes to match input files*/
  /* conattr_all.txt: conattrA conattrP conattrAP */
  const uint_t conattrA_index = 0;
  const uint_t conattrP_index = 1;
  const uint_t conattrAP_index = 2;


#define NUM_FUNCS 8
  uint_t n_total = NUM_FUNCS, n_attr = NUM_FUNCS;
  uint_t attr_indices[NUM_FUNCS];
  static double lambda_values[NUM_FUNCS]; /* init to zero, unused */
  double exponent_values[NUM_FUNCS], obs_stats[NUM_FUNCS];
  static double theta[NUM_FUNCS]; /* init to zero, unused */
  attr_change_stats_func_t *attr_change_stats_funcs[NUM_FUNCS];


  /*
   * First do statisics involving categorical attributes, outputting
   * statistic values for comparison against statnet output from
   * b1nodematch and b2nodematch (see script
   * run_test_b1nodematch_bpnet_A12000_B4000_attr.sh)
   */
  
  attr_change_stats_funcs[0] = &changeBipartiteNodematchAlphaA;
  attr_indices[0]            = catattrA_index;
  exponent_values[0]         = exponent;
  
  attr_change_stats_funcs[1] = &changeBipartiteNodematchBetaA;
  attr_indices[1]            = catattrA_index;
  exponent_values[1]         = exponent;

  attr_change_stats_funcs[2] = &changeBipartiteNodematchAlphaB;
  attr_indices[2]            = catattrP_index;
  exponent_values[2]         = exponent;
  
  attr_change_stats_funcs[3] = &changeBipartiteNodematchBetaB;
  attr_indices[3]            = catattrP_index;
  exponent_values[3]         = exponent;
  
  attr_change_stats_funcs[4] = &changeBipartiteNodematchAlphaA;
  attr_indices[4]            = catattrAP_index;
  exponent_values[4]         = exponent;
  
  attr_change_stats_funcs[5] = &changeBipartiteNodematchBetaA;
  attr_indices[5]            = catattrAP_index;
  exponent_values[5]         = exponent;

  attr_change_stats_funcs[6] = &changeBipartiteNodematchAlphaB;
  attr_indices[6]            = catattrAP_index;
  exponent_values[6]         = exponent;
  
  attr_change_stats_funcs[7] = &changeBipartiteNodematchBetaB;
  attr_indices[7]            = catattrAP_index;
  exponent_values[7]         = exponent;

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
  assert(DOUBLE_APPROX_EQ(stat_value,  obs_stats[0]));

  stat_value= BipartiteNodematchBetaA(g, attr_indices[1], exponent_values[1]);
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[1]));

  stat_value= BipartiteNodematchAlphaB(g, attr_indices[2], exponent_values[2]);
  assert(DOUBLE_APPROX_EQ(stat_value,  obs_stats[2]));

  stat_value= BipartiteNodematchBetaB(g, attr_indices[3], exponent_values[3]);
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[3]));

  stat_value= BipartiteNodematchAlphaA(g, attr_indices[4], exponent_values[4]);
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[4]));

  stat_value= BipartiteNodematchBetaA(g, attr_indices[5], exponent_values[5]);
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[5]));

  stat_value= BipartiteNodematchAlphaB(g, attr_indices[6], exponent_values[6]);
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[6]));

  stat_value= BipartiteNodematchBetaB(g, attr_indices[7], exponent_values[7]);
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[7]));

  free_graph(g);

#undef NOT_USED
#ifdef NOT_USED
  /* experimntal statisic BipartiteDiffBetaA not used as no
     change statistic implemented */
  
  /*
   * Now do statistics involving continuous attribute, but no output
   * from this as there is no equivalent in statnet to test against,
   * we just use the internal tests (assertions) here.
   */
  g = allocate_graph(num_nodes, FALSE/*undirected*/, TRUE/*bipartite*/,
		     num_A_nodes);

  if (load_attributes(g, NULL, catattr_filename, conattr_filename,
		      NULL) != 0) {
    fprintf(stderr, "ERRROR: load node attributes failed\n");
    exit(1);
  }

  if (!(file = fopen(edgelist_filename, "r"))) {
    fprintf(stderr, "error opening file %s (%s)\n",
            edgelist_filename, strerror(errno));
    return -1;
  }


  n_total = 1;
  n_attr = 1;

  attr_change_stats_funcs[0] = &changeBipartiteDiffBetaA;
  attr_indices[0]            = conattrA_index;
  exponent_values[0]         = exponent;

  for (i = 0; i < n_total; i++) {
    obs_stats[i] = 0;
  }
  g = load_graph_from_arclist_file(file, g, TRUE,
                                   n_total, n_attr, 0, 0, NULL,
                                   lambda_values, attr_change_stats_funcs,
                                   NULL, NULL, attr_indices, exponent_values,
                                   NULL, obs_stats, theta);

  stat_value= BipartiteDiffBetaA(g, attr_indices[0], exponent_values[0]);
  assert(DOUBLE_APPROX_EQ(stat_value,  obs_stats[0]));

  free_graph(g);
  
#endif /*NOT_USED*/
  
  exit(0);
}
