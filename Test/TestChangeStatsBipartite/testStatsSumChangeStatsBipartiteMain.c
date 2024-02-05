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
 * Usage:  testStatsSumChangeStatsBipartite  [-s]<in_edgelistfile> <lambda>
 *          -s : also test with slow implementations of statistic funcion
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
#include <getopt.h>
#include "graph.h"
#include "changeStatisticsUndirected.h"
#include "changeStatisticsBipartiteUndirected.h"
#include "loadGraph.h"



/*****************************************************************************
 *
 * utility functions
 *
 ****************************************************************************/

/* Approximate double floating point equality */
#define DOUBLE_APPROX_EQ_TEST(a, b) ( fabs((a) - (b)) <= 1e-06 )

/*
 * Binomial coefficient n choose 2
 */
static ulong_t n_choose_2(uint_t n)
{
  if (n < 2) {
    return 0;
  }
  return n * (n - 1) / 2;
}

/*
 * Binomial coefficient n choose k
 */
static double n_choose_k(uint_t n, uint_t k)
{
  uint_t i;
  double a = 1, b = 1;
  uint_t l = k;

  if (n < k) {
    return 0;
  }

  if (n - k < k) {
    l = n - k;
  }
  
  for (i = 1; i <= l; i++) {
    a *= (n + 1 - i);
    b *= i;
  }
  /*printf("%u %u %llu %llu %llu\n", n, k, a, b, a/b);*/
  //assert(a % b == 0);
  return a / b;
}

/*
 * count k-two-paths, as defined by eqn (6.11) in
 *
 *   Wang, P., Sharpe, K., Robins, G. L., & Pattison,
 *   P. E. (2009). Exponential random graph (p∗) models for affiliation
 *   networks. Social Networks, 31(1), 12-25.
 *
 */
static ulonglong_t k_two_paths_A(const graph_t *g, uint_t k)
{
  uint_t l,i;
  ulonglong_t count = 0;

  assert(k > 0);

  for (l = g->num_A_nodes + 1; l < g->num_A_nodes + g->num_B_nodes; l++){
    for (i = g->num_A_nodes; i < l; i++) {
      assert(bipartite_node_mode(g, i) == MODE_B);
      assert(bipartite_node_mode(g, l) == MODE_B);
      count += n_choose_k(GET_B2PATH_ENTRY(g, i, l), k);
    }
  }
  /* Note despite eqn (6.11) in Wang et al. (2009) having this
     statistic multiplied by 1/2 when k = 2 "due to symmetry", this is
     not actually correct for bipartite networks here as we are
     considering only either mode A or mdoe B nodes (not both at
     once). So the symmetry which exists for one-mode networks as per
     Snijders et al. (2006) "New specifications for exponential random
     graph models" [eqns (25a,b) and (26a), pp. 123-124] when the
     summation over all i < j means two pairs of nodes are considered
     in a four-cycle, is not true here as the summation over i < l
     considers only one pair of nodes (mode B nodes in k_two_paths_A()
     or mode A nodes in k_two_paths_B) are considered.
  */
  return count;
}

static ulonglong_t k_two_paths_B(const graph_t *g, uint_t k)
{
  uint_t l,i;
  ulonglong_t count = 0;

  assert(k > 0);

  for (l = 1; l < g->num_A_nodes; l++){
    for (i = 0; i < l; i++) {
      assert(bipartite_node_mode(g, i) == MODE_A);
      assert(bipartite_node_mode(g, l) == MODE_A);
      count += n_choose_k(GET_A2PATH_ENTRY(g, i, l), k);
    }
  }
  /* Note despite eqn (6.11) in Wang et al. (2009) having this
     statistic multiplied by 1/2 when k = 2 "due to symmetry", this is
     not actually correct for bipartite networks here as we are
     considering only either mode A or mdoe B nodes (not both at
     once). So the symmetry which exists for one-mode networks as per
     Snijders et al. (2006) "New specifications for exponential random
     graph models" [eqns (25a,b) and (26a), pp. 123-124] when the
     summation over all i < j means two pairs of nodes are considered
     in a four-cycle, is not true here as the summation over i < l
     considers only one pair of nodes (mode B nodes in k_two_paths_A()
     or mode A nodes in k_two_paths_B) are considered.
  */
  return count;
}


/*****************************************************************************
 *
 * statistics functions (summing change statistics is verified against these)
 *
 ****************************************************************************/

/*
 * Statistic for FourCycles, number of four-cycles in a bipartite graph.
 *
 * This version counting over pairs of mode A nodes, but result must be
 * equal to that counting over pairs of mode B nodes instead.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      number of four-cycles in bipartite graph g
 */
static double FourCyclesA(const graph_t *g)
{
  uint_t i,l;
  double value = 0;

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = 1; i < g->num_A_nodes; i++) {
    for (l = 0; l < i; l++) {
      assert(bipartite_node_mode(g, i) == MODE_A);
      assert(bipartite_node_mode(g, l) == MODE_A);
      value += n_choose_k(GET_A2PATH_ENTRY(g, i, l), 2);
    }
  }
  return value;
}

/*
 * Statistic for FourCycles, number of four-cycles in a bipartite graph.
 *
 * This version counting over pairs of mode B nodes, but result must be
 * equal to that counting over pairs of mode A nodes instead.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      number of four-cycles in bipartite graph g
 */
static double FourCyclesB(const graph_t *g)
{
  uint_t i,l;
  double value = 0;

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = g->num_A_nodes + 1; i < g->num_A_nodes + g->num_B_nodes; i++) {
    for (l = g->num_A_nodes; l < i; l++) {
      assert(bipartite_node_mode(g, i) == MODE_B);
      assert(bipartite_node_mode(g, l) == MODE_B);
      value += n_choose_k(GET_B2PATH_ENTRY(g, i, l), 2);
    }
  }
  return value;
}


/*
 * Statistic for BipartiteAltKCyclesA, alternating k-cycles for type A
 * nodes (K-Ca in BPNet, XACA in MPNet) defined by eqn (6.12) in:
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
static double BipartiteAltKCyclesA(const graph_t *g, double lambda)
{
  uint_t i,l;
  double value = 0;

  assert (lambda > 1.0);

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (l = g->num_A_nodes + 1; l < g->num_A_nodes + g->num_B_nodes; l++) {
    for (i = g->num_A_nodes; i < l; i++) {
      assert(bipartite_node_mode(g, i) == MODE_B);
      assert(bipartite_node_mode(g, l) == MODE_B);
      value += 1 - POW_LOOKUP(1-1/lambda, GET_B2PATH_ENTRY(g, i, l));
    }
  }
  return lambda * value;
}


/*
 * Statistic for BipartiteAltKCyclesB, alternating k-cycles for type B
 * nodes (K-Cp in BPNet, XACB in MPNet) defined by eqn (6.12) in:
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
static double BipartiteAltKCyclesB(const graph_t *g, double lambda)
{
  uint_t i,l;
  double value = 0;

  assert (lambda > 1.0);

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (l = 1; l < g->num_A_nodes; l++) {
    for (i = 0; i < l; i++) {
      assert(bipartite_node_mode(g, i) == MODE_A);
      assert(bipartite_node_mode(g, l) == MODE_A);
      value += 1 - POW_LOOKUP(1-1/lambda, GET_A2PATH_ENTRY(g, i, l));
    }
  }
  return lambda * value;
}




/*
 * Statistic for BipartiteAltKCyclesA, alternating k-cycles for type A
 * nodes (K-Ca in BPNet, XACA in MPNet), alternative (inefficient)
 * implementation, defined by eqn (6.12) in:
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
static double BipartiteAltKCyclesA_SLOW(const graph_t *g, double lambda)
{
  uint_t i;
  double value;

  assert (lambda > 1.0);

  assert(g->is_bipartite);
  assert(!g->is_directed);

  /* Note despite eqn (6.12) in Wang et al. (2009) multipliying the
     second term [k_two_paths_A(g, 2)/lambda] by 2, this is not
     actually correct as the division by two in eqn (6.11) is not
     correct (see comment in k_two_paths_A()), so there is no factor
     of 2 here.
  */
  value = k_two_paths_A(g, 1) - k_two_paths_A(g, 2)/lambda;

  for (i = 3; i < g->num_A_nodes + g->num_B_nodes - 1; i++) {
    value += pow(-1/lambda, i-1) * k_two_paths_A(g, i);
  }
  return value;
}

/*
 * Statistic for BipartiteAltKCyclesB, alternating k-cycles for type B
 * nodes (K-Cp in BPNet, XACB in MPNet), alternative (inefficient)
 * implementation, defined by eqn (6.12) in:
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
static double BipartiteAltKCyclesB_SLOW(const graph_t *g, double lambda)
{
  uint_t i;
  double value;

  assert (lambda > 1.0);

  assert(g->is_bipartite);
  assert(!g->is_directed);

  /* Note despite eqn (6.12) in Wang et al. (2009) multipliying the
     second term [k_two_paths_A(g, 2)/lambda] by 2, this is not
     actually correct as the division by two in eqn (6.11) is not
     correct (see comment in k_two_paths_B()), so there is no factor
     of 2 here.
  */
    value = k_two_paths_B(g, 1) - k_two_paths_B(g, 2)/lambda;

  for (i = 3; i < g->num_A_nodes + g->num_B_nodes - 1; i++) {
    value += pow(-1/lambda, i-1) * k_two_paths_B(g, i);
  }
  return value;
}



/*
 * Statistic for alternating k-4-cycles for type A nodes (new change
 * statistic suggested in email (basically paper outline, with
 * spreadsheet attachments for literature search, examples, etc.)
 * "idea for a (slightly) new bipartite change statistic" sent 23 Nov
 * 2022):
 *
 *   The proposed new statistic is a very simple modification of the
 *   "Alternating k-two-paths" bipartite statistic (K-Ca and K-Cp)
 *   statistics described by Wang et al. (2009, p.19). I propose to
 *   simply remove the first term of Wang et al. (2009) equation 6.12,
 *   and reverse the signs, so that it no longer counts open two-paths,
 *   but the first, positive, term actually counts four-cycles.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *     lambda - decay value > 1.0
 *
 * Return value:
 *      statistic for g with decay value lambda
 */
static double BipartiteAltK4CyclesA_SLOW(const graph_t *g, double lambda)
{
  uint_t i;
  double value;

  assert (lambda > 1.0);

  assert(g->is_bipartite);
  assert(!g->is_directed);

  value = k_two_paths_A(g, 2)/lambda;

  for (i = 3; i < g->num_A_nodes + g->num_B_nodes - 1; i++) {
    value += -1 * pow(-1/lambda, i-1) * k_two_paths_A(g, i);
  }
  return value;
}


/*
 * Count number of four-cycles that a particular node u is involved in.
 *
 * This version for bipartite networks.
 */
static uint_t num_four_cycles_node_SLOW(const graph_t *g, uint_t u)
{
  uint_t v;
  uint_t count = 0;

  if (g->is_bipartite) {
    if (bipartite_node_mode(g, u) == MODE_A) {
      for (v = 0; v < g->num_A_nodes; v++) {
        if (v != u) {
          assert(bipartite_node_mode(g, v) == MODE_A);
          count += n_choose_2(GET_A2PATH_ENTRY(g, u, v));
        }
      }
    } else {
      for (v = g->num_A_nodes; v < g->num_A_nodes + g->num_B_nodes; v++) {
        if (v != u) {
          assert(bipartite_node_mode(g, v) == MODE_B);
          count += n_choose_2(GET_B2PATH_ENTRY(g, u, v));
        }
      }
    }
  }
  return count;
}


/*
 * Statistic for 4-cycles raised to a power. The lambda parameter (>
 * 1.0) (mis)used to specify the value 1/lambda as the epxonent. Note
 * this is not the same meaning of lambda as its original use in the
 * "alternating" parameters.
 *
 * This version counting over pairs of mode A nodes only.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      
 */
static double PowerFourCyclesA(const graph_t *g, double lambda)
{
  uint_t  i;
  double  alpha = 1/lambda;
  double  value = 0;
  uint_t  fourcycle_count = 0;

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = 0; i < g->num_A_nodes; i++) {
    uint_t fourcycle_count_SLOW = num_four_cycles_node_SLOW(g, i);
    fourcycle_count = num_four_cycles_node(g, i);
    assert(fourcycle_count == fourcycle_count_SLOW);
    value += pow(num_four_cycles_node(g, i), alpha);
  }
  return value;
}


/*
 * Statistic for 4-cycles raised to a power. The lambda parameter (>
 * 1.0) (mis)used to specify the value 1/lambda as the epxonent. Note
 * this is not the same meaning of lambda as its original use in the
 * "alternating" parameters.
 *
 * This version counting over pairs of mode B nodes only.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      
 */
static double PowerFourCyclesB(const graph_t *g, double lambda)
{
  uint_t  i;
  double  alpha = 1/lambda;
  double  value = 0;
  uint_t  fourcycle_count = 0;

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = g->num_A_nodes; i < g->num_A_nodes + g->num_B_nodes; i++) {
    uint_t fourcycle_count_SLOW = num_four_cycles_node_SLOW(g, i);
    fourcycle_count = num_four_cycles_node(g, i);
    assert(fourcycle_count == fourcycle_count_SLOW);
    value += pow(num_four_cycles_node(g, i), alpha);
  }
  return value;
}



/*
 * Statistic for 4-cycles raised to a power. The lambda parameter (>
 * 1.0) (mis)used to specify the value 1/lambda as the epxonent. Note
 * this is not the same meaning of lambda as its original use in the
 * "alternating" parameters.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      number of four-cycles in bipartite graph g
 */
static double PowerFourCycles(const graph_t *g, double lambda)
{
  assert(g->is_bipartite);
  assert(!g->is_directed);

  double power_4cycles_A = PowerFourCyclesA(g, lambda);
  double power_4cycles_B = PowerFourCyclesB(g, lambda);
  return power_4cycles_A + power_4cycles_B;
}



/*****************************************************************************
 *
 * main
 *
 ****************************************************************************/

static void usage(const char *progname)
{
  fprintf(stderr, "Usage: %s [-s] <edgelist_file> <lambda>\n"
          "  -s : also test with slow statistics functions\n"
          , progname);
  exit(1);
}

int main(int argc, char *argv[])
{
  int c;
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
  bool also_use_slow_functions = FALSE;
 
  srand(time(NULL));

  while ((c = getopt(argc, argv, "s")) != -1)  {
    switch (c)   {
      case 's':
        also_use_slow_functions = TRUE;
        break;
      default:
        usage(argv[0]);
        break;
    }
  }

  if (argc - optind != 2) {
    usage(argv[0]);
  }
  
  edgelist_filename = argv[optind];
  lambda = strtod(argv[optind+1], &endptr);
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



#define NUM_FUNCS 5
  uint_t n_total = NUM_FUNCS;
  static double lambda_values[NUM_FUNCS];
  double obs_stats[NUM_FUNCS];
  static double theta[NUM_FUNCS]; /* init to zero, unused */
  change_stats_func_t *change_stats_funcs[NUM_FUNCS];


  change_stats_funcs[0] = &changeBipartiteAltKCyclesA;
  lambda_values[0]      = lambda;

  change_stats_funcs[1] = &changeBipartiteAltKCyclesB;
  lambda_values[1]      = lambda;

  change_stats_funcs[2] = &changeFourCycles;
  lambda_values[2]      = 0; /* not used */

  change_stats_funcs[3] = &changePowerFourCycles;
  lambda_values[3]      = lambda;

  change_stats_funcs[4] = &changeBipartiteAltK4CyclesA;
  lambda_values[4]      = lambda;

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

  stat_value = BipartiteAltKCyclesA(g, lambda_values[0]);
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[0]));
  if (also_use_slow_functions) {
    stat_value = BipartiteAltKCyclesA_SLOW(g, lambda_values[0]);
    /* fprintf(stderr,"stat_value   = %.10f\nobs_stats[0] = %.10f\n", stat_value, obs_stats[0]); */
    /* fprintf(stderr, "diff = %g\n", fabs((stat_value) - (obs_stats[0]))); */
    assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[0]));
  }
  

  stat_value = BipartiteAltKCyclesB(g, lambda_values[1]);
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[1]));
  if (also_use_slow_functions) {
    stat_value = BipartiteAltKCyclesB_SLOW(g, lambda_values[1]);
    /* fprintf(stderr,"stat_value   = %.10f\nobs_stats[1] = %.10f\n", stat_value, obs_stats[1]); */
    /* fprintf(stderr, "diff = %g\n", fabs((stat_value) - (obs_stats[1]))); */
    assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[1]));
  }

  stat_value = FourCyclesA(g);
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[2]));
  stat_value = FourCyclesB(g);
  /* fprintf(stderr,"stat_value   = %.10f\nobs_stats[2] = %.10f\n", stat_value, obs_stats[2]); */
  /* fprintf(stderr, "diff = %g\n", fabs((stat_value) - (obs_stats[2]))); */
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[2]));


  stat_value = PowerFourCycles(g, lambda_values[3]);
  fprintf(stderr,"stat_value   = %.10f\nobs_stats[3] = %.10f\n", stat_value, obs_stats[3]);
  fprintf(stderr, "diff = %g\n", fabs((stat_value) - (obs_stats[3])));
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[3]));


  /* test for experimental BipartiteAltK4CyclesA statisic is disabled as the change statisic is not correct */
  const bool DISABLED_TEST_CASE = TRUE;
  if (!DISABLED_TEST_CASE && also_use_slow_functions) {
    stat_value = BipartiteAltK4CyclesA_SLOW(g, lambda_values[4]);
    fprintf(stderr,"stat_value   = %.10f\nobs_stats[4] = %.10f\n", stat_value, obs_stats[4]);
    fprintf(stderr, "diff = %g\n", fabs((stat_value) - (obs_stats[4])));
    assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[4]));
  }

  free_graph(g);
  exit(0);
}
