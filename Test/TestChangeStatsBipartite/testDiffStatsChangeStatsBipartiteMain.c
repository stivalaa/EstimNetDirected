/*****************************************************************************
 * 
 * File:    testDiffStatsChangeStatsBipartiteMain.c
 * Author:  Alex Stivala
 * Created: June 2024
 *
 * Test change statistics implementations by comparing change stats to
 * difference in statistic value computed according to the definition
 * of the statistic (implemented in bipartiteStats.c).
 *
 *
 * Usage:  testDiffStatsChangeStatsBipartite [-s] <in_edgelistfile> <lambda> [nodenums]
 *          -s : also test with slow implementations of statistic funcion
 *
 * Reads graph from Pajek format <in_edgelistfile> and compute stats with
 * weighting parameter <lambda> where real value lambda > 1.
 * If optional nodenums filename specified, then reads (two-column whitespace
 * separated, one pair per line) pairs of nodes to use from there,
 * otherwise randomly generated.
 *
 * Outputs observed statistics value for the statistics, which are computed
 * by summing the change stats over all edges in the data, and verifies
 * that these indeed sum to the statistic value computed directly (in this
 * code).
 *
 * Example:
 * ./testDiffStatsChangeStatsBipartite ../../examples/bipartite/simulated/bpnet_A6000_B750_sparse_sim100000000.net 2.0
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

#define DEFAULT_NUM_TESTS 100


static void usage(const char *progname)
{
  fprintf(stderr, "Usage: %s [-s] <edgelist_file> <lambda> [nodenums]\n"
          "  -s : also test with slow statistics functions\n"
          , progname);
  exit(1);
}



int main(int argc, char *argv[]) 
{
  char buf[1024];
  uint_t i,j;
  char *edgelist_filename = NULL;
  FILE *file           = NULL;
  uint_t     num_nodes = 0;
  uint_t     num_A_nodes = 0;
  /*  uint_t     num_P_nodes = 0; */
  graph_t *g         = NULL;
  int    num_tests     = 0;
  int    readNodeNums  = FALSE;
  char  *nodenumfilename = NULL;
  FILE  *nodenumfile     = NULL;
  double lambda;
  char  *endptr; /* for strtod() */
  bool also_use_slow_functions = FALSE;
  int  c;
  bool edge_removed = FALSE;
  double delta_BipartitePowerFourCyclesA, delta_BipartitePowerFourCyclesB,
    delta_PowerFourCycles, delta_BipartiteAltKCyclesA,
    delta_BipartiteAltKCyclesB;
  double without_BipartitePowerFourCyclesA, without_BipartitePowerFourCyclesB,
    without_PowerFourCycles, without_BipartiteAltKCyclesA,
    without_BipartiteAltKCyclesB;


 
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

  if (argc - optind < 2 || argc - optind > 3) {
    usage(argv[0]);
  }
  
  edgelist_filename = argv[optind];
  lambda = strtod(argv[optind+1], &endptr);
  if (lambda <= 1.0) {
    fprintf(stderr, "lambda value %g is not > 1.0\n", lambda);
    return -1;
  }
  if (argc - optind == 3) {
    readNodeNums = TRUE;
    nodenumfilename = argv[optind+2];
    if (!(nodenumfile = fopen(nodenumfilename, "r"))) {
      fprintf(stderr, "open %s for read failed (%s)\n", nodenumfilename,
              strerror(errno));
      exit(1);
    }
  } else {
    readNodeNums = FALSE;
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
  g = load_graph_from_arclist_file(file, g, FALSE,
                                     0, 0, 0, 0, NULL, NULL, NULL, NULL,
                                     NULL, NULL, NULL, NULL, NULL, NULL);

  num_tests = 0;
  while (TRUE) {
    if (readNodeNums) {
      if (feof(nodenumfile)) {
        break;
      }
      if (!fgets(buf, sizeof(buf)-1, nodenumfile)) {
        break;
      }
      if (sscanf(buf, "%d %d\n", &i, &j) != 2) {
        fprintf(stderr, "error reading node nums\n");
        exit(1);
      }
      if (i >= g->num_A_nodes) {
        fprintf(stderr, "bad i node num %d\n", i);
        exit(1);
      }
      if (j >= g->num_B_nodes) {
        fprintf(stderr, "bad j node num %d\n", j);
        exit(1);
      }
      j += g->num_A_nodes; /* in file is 0 .. N-Np-1 but we number Np .. N-1 */
    } else {
      i = rand() % g->num_A_nodes;
      j = g->num_A_nodes + rand() % g->num_B_nodes;
    }

    if (isEdge(g, i, j)) {
      removeEdge(g, i, j);
      edge_removed = TRUE;
    }

    delta_BipartitePowerFourCyclesA = changeBipartitePowerFourCyclesA(g, i, j,
                                                                      lambda);
    delta_BipartitePowerFourCyclesB = changeBipartitePowerFourCyclesB(g, i, j,
                                                                      lambda);
    delta_PowerFourCycles = changePowerFourCycles(g, i, j, lambda);
    delta_BipartiteAltKCyclesA = changeBipartiteAltKCyclesA(g, i, j, lambda);
    delta_BipartiteAltKCyclesB = changeBipartiteAltKCyclesB(g, i, j, lambda);

    /* verify that the sum of two-mode changBipartitePowerFourCyclesA and
       changeBipartiteFourCyclesB is equal to the one-mode
       changePowerFourCycles applied to a bipartite graph */
    assert(DOUBLE_APPROX_EQ_TEST(delta_BipartitePowerFourCyclesA +
                                 delta_BipartitePowerFourCyclesB,
                                 delta_PowerFourCycles));

    /* verify that change statisic is used to difference of statistic
       computed with edge and without edge */
    without_BipartitePowerFourCyclesA = PowerFourCyclesA(g, lambda);
    without_BipartitePowerFourCyclesB = PowerFourCyclesB(g, lambda);
    without_BipartiteAltKCyclesA = BipartiteAltKCyclesA(g, lambda);
    if (also_use_slow_functions) {
      assert(DOUBLE_APPROX_EQ_TEST(BipartiteAltKCyclesA_SLOW(g, lambda),
                                   without_BipartiteAltKCyclesA));
    }
    insertEdge(g, i, j);
    assert(DOUBLE_APPROX_EQ_TEST(delta_BipartitePowerFourCyclesA,
                                 PowerFourCyclesA(g, lambda) -
                                 without_BipartitePowerFourCyclesA));
    assert(DOUBLE_APPROX_EQ_TEST(delta_BipartitePowerFourCyclesB,
                                 PowerFourCyclesB(g, lambda) -
                                 without_BipartitePowerFourCyclesB));
    assert(DOUBLE_APPROX_EQ_TEST(delta_BipartiteAltKCyclesA,
                                 BipartiteAltKCyclesA(g, lambda) -
                                 without_BipartiteAltKCyclesA));
    removeEdge(g, i, j);

    num_tests++;
    if (edge_removed) {
      insertEdge(g, i, j);
    }
    if (!readNodeNums && num_tests >= DEFAULT_NUM_TESTS) {
      break;
    }
  }
  if (readNodeNums) {
    fclose(nodenumfile);
    if (!(nodenumfile = fopen(nodenumfilename, "r"))) {
      fprintf(stderr, "open %s for read failed (%s)\n", nodenumfilename,
              strerror(errno));
      exit(1);
    }
  }
  
  
  if (readNodeNums) {
    fclose(nodenumfile);
  }
  free_graph(g);
  exit(0);
}
