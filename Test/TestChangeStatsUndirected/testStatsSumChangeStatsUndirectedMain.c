/*****************************************************************************
 * 
 * File:    testStatsSumChangeStatsUndirectedMain.c
 * Author:  Alex Stivala
 * Created: January 2024
 *
 * Test change statistics implementations by comparing sum of change
 * stats for all edges in network to statistic value computed according
 * to the definition of the statistic (implemented in undirectedStats.c for
 * this test module).
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
#include "undirectedStats.h"



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
  ulonglong_t stat_value_int, alt_stat_value_int;
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

  stat_value_int = FourCycles(g);
  alt_stat_value_int = FourCycles_sum_by_node(g);
  fprintf(stderr, "stat_value_int = %llu\n", stat_value_int);/*XXX*/
  assert(stat_value_int ==   alt_stat_value_int);
  assert(DOUBLE_APPROX_EQ_TEST((double)stat_value_int,  obs_stats[0]));

  stat_value = PowerFourCycles(g, lambda_values[1]);
  fprintf(stderr,"stat_value   = %.10f\nobs_stats[1] = %.10f\n", stat_value, obs_stats[1]);
  fprintf(stderr, "diff = %g\n", fabs((stat_value) - (obs_stats[1])));
  assert(DOUBLE_APPROX_EQ_TEST(stat_value,  obs_stats[1]));

  

  free_graph(g);
  exit(0);
}
