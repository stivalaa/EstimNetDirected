/*****************************************************************************
 * 
 * File:    testBoundaryCrossingRatio.c
 * Author:  Alex Stivala
 * Created: March 2020
 *
 * Test boundary crossing ratio computation.
 *
 *
 * Usage:  testBoundaryCrossingRatio arclistFile setattrFile
 *
 * Reads digraph in Pajek format from arclistFile and set attributes from
 * setattrFile. Write output to stdout in format:
 *
 * node ratio
 *
 * where node is the node number and ratio is the boundary crossing ratio
 * value for that node.
 *
 ****************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "digraph.h"
#include "changeStatisticsDirected.h"
#include "loadDigraph.h"

int main(int argc, char *argv[]) 
{
  uint_t i;
  char *arclist_filename = NULL;
  char *set_filename = NULL;
  FILE *file           = NULL;
  uint_t     num_nodes = 0;
  digraph_t *g         = NULL;
  struct timeval start_timeval, end_timeval, elapsed_timeval;
  int    etime;
  

  srand(time(NULL));

  if (argc != 3) {
    fprintf(stderr, "Usage: %s <arclist_filename> <setattr_filename>\n",
            argv[0]);
    exit(1);
  }
  arclist_filename = argv[1];
  set_filename = argv[2];
  if (!(file = fopen(arclist_filename, "r"))) {
    fprintf(stderr, "error opening file %s (%s)\n", 
            arclist_filename, strerror(errno));
    exit(1);
  }
  gettimeofday(&start_timeval, NULL);
  fprintf(stderr, "loading arc list from %s and building two-path tables...",
         arclist_filename);
  num_nodes = get_num_vertices_from_arclist_file(file); /* closes file */
  g = allocate_digraph(num_nodes);
  if (!(file = fopen(arclist_filename, "r"))) {
    fprintf(stderr, "error opening file %s (%s)\n", 
            arclist_filename, strerror(errno));
    exit(1);
  }
  g = load_digraph_from_arclist_file(file, g, FALSE,
                                     0, 0, 0, 0, NULL, NULL, NULL, NULL,
                                     NULL, NULL, NULL, NULL);
  gettimeofday(&end_timeval, NULL);
  timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
  etime = 1000 * elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
  fprintf(stderr, "%.2f s\n", (double)etime/1000);
#ifdef DEBUG_DIGRAPH
  dump_digraph_arclist(g);
#endif /*DEBUG_DIGRAPH*/
  
  if (load_attributes(g, NULL, NULL, NULL, set_filename)) {
    fprintf(stderr, "error loading set attributes file %s\n", set_filename);
    exit(1);
  }
  
  for (i = 0; i < g->num_nodes; i++) {
    printf("%u %g\n", i, boundary_crossing_ratio(g, i, 0));
  }
  exit(0);
}
  
