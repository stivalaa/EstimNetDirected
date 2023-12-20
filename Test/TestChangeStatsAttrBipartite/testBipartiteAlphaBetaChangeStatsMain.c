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
  char buf[1024];
  uint_t i,j;
  char *edgelist_filename = NULL;
  FILE *file           = NULL;
  uint_t     num_nodes = 0;
  uint_t     num_A_nodes = 0;
  /*  uint_t     num_P_nodes = 0; */
  graph_t *g         = NULL;
  struct timeval start_timeval, end_timeval, elapsed_timeval;
  int    etime;
  char *binattr_filename, *conattr_filename, *catattr_filename;
 
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
  
  gettimeofday(&start_timeval, NULL);
#ifdef TWOPATH_LOOKUP
  fprintf(stderr, "loading edge list from %s and building two-path tables...",
         edgelist_filename);
#else
    fprintf(stderr, "loading edge list from %s...",
         edgelist_filename);
#endif
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
  g = load_graph_from_arclist_file(file, g, FALSE,
                                     0, 0, 0, 0, NULL, NULL, NULL, NULL,
                                     NULL, NULL, NULL, NULL, NULL, NULL);
  gettimeofday(&end_timeval, NULL);
  timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
  etime = 1000 * elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
  fprintf(stderr, "%.2f s\n", (double)etime/1000);

  
  /* hardcoding indices of attributes to match input files*/
  /* catattr_all.txt: catattrA catattrP catattrAP */
  uint_t catattrA_index = 0;
  uint_t catattrP_index = 1;
  uint_t catattrAP_index = 2;
  
    
  free_graph(g);
  exit(0);
}
