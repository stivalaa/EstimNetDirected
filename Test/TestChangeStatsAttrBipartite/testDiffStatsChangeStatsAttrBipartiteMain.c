/*****************************************************************************
 * 
 * File:    testDiffStatsChangeStatsAttrBipartiteMain.c
 * Author:  Alex Stivala
 * Created: January 2025
 *
 * Test change statistics implementations by comparing change stats to
 * difference in statistic value computed according to the definition
 * of the statistic (implemented in attrbipartiteStats.c).
 *
 *
 * Usage:  testDiffStatsChangeStatsAttrBipartite <in_edgelistfile> <binattrfile> [nodenums]
 *
 * Reads graph from Pajek format <in_edgelistfile> and binary
 * attributes from <binattrfile>.  If optional nodenums filename
 * specified, then reads (two-column whitespace separated, one pair
 * per line) pairs of nodes to use from there, otherwise randomly
 * generated.
 *
 * Example:
 * ./testDiffStatsChangeAttrStatsBipartite  ../../examples/bipartite/simulated/bpnet_A12000_B4000_attrs_sim830000000.net  ../../examples/bipartite/simulation/binattr_all.txt
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
#include "attrbipartiteStats.h"

#define DEFAULT_NUM_TESTS 100


static void usage(const char *progname)
{
  fprintf(stderr, "Usage: %s <edgelist_file> <binattr_file> [nodenums]\n",
          progname);
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
  bool edge_removed = FALSE;
  char *binattr_filename;
  double delta_BipartiteExactlyOneNeighbourA,
    delta_BipartiteExactlyOneNeighbourB;
  double without_BipartiteExactlyOneNeighbourA,
    without_BipartiteExactlyOneNeighbourB;
 
  srand(time(NULL));

  if (argc  < 3 || argc  > 4) {
    usage(argv[0]);
  }
  
  edgelist_filename = argv[1];
  binattr_filename = argv[2];
  if (argc  == 4) {
    readNodeNums = TRUE;
    nodenumfilename = argv[3];
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
  if (load_attributes(g, binattr_filename, NULL, NULL, NULL) != 0) {
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


  
  /* hardcoding indices of attributes to match input files*/
  /* binattr_all.txt: binattrA binattrP binattrAP */
  const uint_t binattrA_index = 0;
  const uint_t binattrP_index = 1;
  const uint_t binattrAP_index = 2;

  

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

    delta_BipartiteExactlyOneNeighbourA = changeBipartiteExactlyOneNeighbourA(g, i, j, binattrP_index, FALSE, 0);
    


    /* verify that change statisic is equal to difference of statistic
       computed with edge and without edge */
    without_BipartiteExactlyOneNeighbourA = BipartiteExactlyOneNeighbourA(g, binattrP_index);

    insertEdge(g, i, j);
    
    assert(DOUBLE_APPROX_EQ_TEST(delta_BipartiteExactlyOneNeighbourA,
                                 BipartiteExactlyOneNeighbourA(g, binattrP_index) -
                                 without_BipartiteExactlyOneNeighbourA));
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
  }
  free_graph(g);
  exit(0);
}
