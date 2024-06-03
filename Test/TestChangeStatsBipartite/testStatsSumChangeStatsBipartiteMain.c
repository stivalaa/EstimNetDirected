/*****************************************************************************
 * 
 * File:    testStatsSumChangeStatsBipartiteMain.c
 * Author:  Alex Stivala
 * Created: January 2024
 *
 * Test change statistics implementations by comparing sum of change
 * stats for all edges in network to statistic value computed according
 * to the definition of the statistic (implemented in bipartiteStats.c).
 *
 *
 * Usage:  testStatsSumChangeStatsBipartite [-s] <in_edgelistfile> <lambda>
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
#include "bipartiteStats.h"



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



#define NUM_FUNCS 7
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

  change_stats_funcs[4] = &changeBipartitePowerFourCyclesA;
  lambda_values[4]      = lambda;

  change_stats_funcs[5] = &changeBipartitePowerFourCyclesB;
  lambda_values[5]      = lambda;

  change_stats_funcs[6] = &changeBipartiteAltK4CyclesA;
  lambda_values[6]      = lambda;

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

  stat_value = PowerFourCyclesA(g, lambda_values[4]);
  fprintf(stderr,"stat_value   = %.10f\nobs_stats[4] = %.10f\n", stat_value, obs_stats[4]);
  fprintf(stderr, "diff = %g\n", fabs((stat_value) - (obs_stats[4])));
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[4]));

  stat_value = PowerFourCyclesB(g, lambda_values[5]);
  fprintf(stderr,"stat_value   = %.10f\nobs_stats[5] = %.10f\n", stat_value, obs_stats[5]);
  fprintf(stderr, "diff = %g\n", fabs((stat_value) - (obs_stats[5])));
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[5]));


  /* test for experimental BipartiteAltK4CyclesA statisic is disabled as the change statisic is not correct */
  const bool DISABLED_TEST_CASE = TRUE;
  if (!DISABLED_TEST_CASE && also_use_slow_functions) {
    stat_value = BipartiteAltK4CyclesA_SLOW(g, lambda_values[6]);
    fprintf(stderr,"stat_value   = %.10f\nobs_stats[6] = %.10f\n", stat_value, obs_stats[6]);
    fprintf(stderr, "diff = %g\n", fabs((stat_value) - (obs_stats[6])));
    assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[6]));
  }

  free_graph(g);
  exit(0);
}
