/*****************************************************************************
 * 
 * File:    testStatsSumChangeStatsUndirectedMain.c
 * Author:  Alex Stivala
 * Created: January 2024
 *
 * Test change statistics implementations by comparing sum of change
 * stats for all edges in network to statistic value computed according
 * to the definition of the statistic (implemented in this test module).
 *
 *
 * Usage:  testStatsSumChangeStatsUndirected  [-s]<in_edgelistfile> <lambda>
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
 * ./testStatsSumChangeStatsUndirected netscience_edgelist.txt 2.0
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
static ulonglong_t n_choose_k(uint_t n, uint_t k)
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
  return a / b;
}


/*****************************************************************************
 *
 * statistics functions (summing change statistics is verified against these)
 *
 ****************************************************************************/

/*
 * Statistic for FourCycles, number of four-cycles in an undirected graph.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      number of four-cycles in undirected graph g
 */
static ulonglong_t FourCycles(const graph_t *g)
{
  uint_t i,l;
  ulonglong_t value = 0;

  assert(!g->is_bipartite);
  assert(!g->is_directed);

  for (i = 1; i < g->num_nodes; i++) {
    for (l = 0; l < i; l++) {
      value += n_choose_k(GET_2PATH_ENTRY(g, i, l), 2);
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
static uint_t num_four_cycles_node(const graph_t *g, uint_t u)
{
  uint_t k,v;
  uint_t count = 0;

  /* TODO implement this more efficiently instead of iterating over all nodes */
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
  } else {
    for (v = 0; v < g->num_nodes; v++){
      if (v != u) {
        count += n_choose_2(GET_2PATH_ENTRY(g, u, v));
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
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      
 */
static double PowerFourCycles(const graph_t *g, double lambda)
{
  uint_t  i,l;
  double  alpha = 1/lambda;
  ulong_t count = 0;
  double  value = 0;

  assert(!g->is_directed);

  for (i = 0; i < g->num_nodes; i++) {
    value += pow(num_four_cycles_node(g, i), alpha);
  }
  return value;
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

  (void)also_use_slow_functions; /* unused variable for now */
  
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
  
  num_nodes = get_num_vertices_from_arclist_file(file);/* closes file */
  
  g = allocate_graph(num_nodes, FALSE/*is_directed*/, FALSE/*is_bipartite*/,
                     num_A_nodes);

  if (!(file = fopen(edgelist_filename, "r"))) {
    fprintf(stderr, "error opening file %s (%s)\n",
            edgelist_filename, strerror(errno));
    return -1;
  }



#define NUM_FUNCS 2
  uint_t n_total = NUM_FUNCS;
  static double lambda_values[NUM_FUNCS];
  double obs_stats[NUM_FUNCS];
  static double theta[NUM_FUNCS]; /* init to zero, unused */
  change_stats_func_t *change_stats_funcs[NUM_FUNCS];

  change_stats_funcs[0] = &changeFourCycles;
  lambda_values[0]      = 0; /* not used */

  change_stats_funcs[1] = &changePowerFourCycles;
  lambda_values[1]      = lambda;

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

  stat_value = FourCycles(g);
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[0]));

  stat_value = PowerFourCycles(g, lambda_values[1]);
  fprintf(stderr,"stat_value   = %.10f\nobs_stats[1] = %.10f\n", stat_value, obs_stats[1]);
  fprintf(stderr, "diff = %g\n", fabs((stat_value) - (obs_stats[1])));
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[1]));

  

  free_graph(g);
  exit(0);
}
