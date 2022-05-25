/*****************************************************************************
 * 
 * File:    testChangeStatsBipartite.c
 * Author:  Alex Stivala
 * Created: October 2017
 *
 * Test bipartite change stats and two-path table update.
 *
 *
 * Usage:  testChangeStatsBipartite  <in_edgelistfile> [nodenums]
 *
 * Reads graph from Pajek format <in_edgelistfile>.
 * If optional nodenums filename specified, then reads (two-column whitespace
 * separated, one pair per line) pairs of nodes to use from there,
 * otherwise randomly generated.
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
#include "changeStatisticsUndirected.h"
#include "changeStatisticsBipartiteUndirected.h"
#include "loadGraph.h"

#define DEFAULT_NUM_TESTS 1000

#ifdef TWOPATH_LOOKUP
/* get stats and dump two-path hash table */
static void dumpTwoPathTable(const graph_t *g) {
  uint_t sumA,maxA,nnzA;
  uint_t sumB,maxB,nnzB;
#ifdef TWOPATH_HASHTABLES
  twopath_record_t *r;
  uint_t val;
#else
  uint_t i, j;
#endif /*TWOPATH_HASHTABLES*/
  
  sumA = maxA = sumB = maxB = 0;
  nnzA = nnzB = 0;

#ifdef TWOPATH_HASHTABLES
  nnzA = HASH_COUNT(g->twoPathHashTabA);
  nnzB = HASH_COUNT(g->twoPathHashTabB);

  for (r = g->twoPathHashTabA; r !=NULL; r = r->hh.next) {
    val = r->value;
    sumA += val;
    if (val > maxA) {
      maxA = val;
    }
  }
  for (r = g->twoPathHashTabB; r !=NULL; r = r->hh.next) {
    val = r->value;
    sumB += val;
    if (val > maxB) {
      maxB = val;
    }
  }
#else
 for (i = 0; i < g->num_A_nodes; i++) {
    for (j = 0; j < g->num_A_nodes; j++) {
      sumA += g->twoPathMatrixA[INDEX2D(i, j, g->num_A_nodes)];
      if (g->twoPathMatrixA[INDEX2D(i, j, g->num_A_nodes)] > 0) {
        nnzA++;
      }
      if (g->twoPathMatrixA[INDEX2D(i, j, g->num_A_nodes)] > maxA) {
        maxA = g->twoPathMatrixA[INDEX2D(i, j, g->num_A_nodes)];
      }
    }
  }
 for (i = g->num_A_nodes; i < g->num_A_nodes + g->num_B_nodes; i++) {
    for (j = g->num_A_nodes; j < g->num_A_nodes + g->num_B_nodes; j++) {
      /* Note subtracting num_A_nodes as B nodes are numbered
	 num_A_nodes .. num_nodes */
      sumB += g->twoPathMatrixB[INDEX2D(i-g->num_A_nodes, j-g->num_A_nodes, g->num_B_nodes)];
      if (g->twoPathMatrixB[INDEX2D(i-g->num_A_nodes, j-g->num_A_nodes, g->num_B_nodes)] > 0) {
        nnzB++;
      }
      if (g->twoPathMatrixB[INDEX2D(i-g->num_A_nodes, j-g->num_A_nodes, g->num_B_nodes)] > maxB) {
        maxB = g->twoPathMatrixB[INDEX2D(i-g->num_A_nodes, j-g->num_A_nodes, g->num_B_nodes)];
      }
    }
  }
 #endif /*TWOPATH_HASHTABLES*/

  printf("vP2p sum = %u, max = %u\n", sumA, maxA);
  printf("vA2p sum = %u, max = %u\n", sumB, maxB);
  printf("vP2p nnz = %u (%.4f%%)\n", nnzA,
         100*(double)nnzA/(g->num_A_nodes*g->num_A_nodes));
  printf("vA2p nnz = %u (%.4f%%)\n", nnzB,
         100*(double)nnzB/(g->num_B_nodes*g->num_B_nodes));
}
#endif /*TWOPATH_LOOKUP*/


int main(int argc, char *argv[]) 
{
  const double lambda = 2.0;  /* always use decay lambda = 2.0 on tests */
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
  struct timeval start_timeval, end_timeval, elapsed_timeval;
  int    etime;
 
  srand(time(NULL));

  if (argc < 2 || argc > 3) {
    fprintf(stderr, "Usage: %s <inedgelist_file> [nodenumsfile]\n", argv[0]);
    exit(1);
  }
  edgelist_filename = argv[1];
  if (argc == 3) {
    readNodeNums = TRUE;
    nodenumfilename = argv[2];
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
  if (!(file = fopen(edgelist_filename, "r"))) {
    fprintf(stderr, "error opening file %s (%s)\n", 
            edgelist_filename, strerror(errno));
    return -1;
  }
  g = load_graph_from_arclist_file(file, g, FALSE,
                                     0, 0, 0, 0, NULL, NULL, NULL, NULL,
                                     NULL, NULL, NULL, NULL, NULL);
  gettimeofday(&end_timeval, NULL);
  timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
  etime = 1000 * elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
  fprintf(stderr, "%.2f s\n", (double)etime/1000);
#ifdef DEBUG_DIGRAPH
  dump_graph_arclist(g);
#endif /*DEBUG_DIGRAPH*/

#ifdef TWOPATH_LOOKUP
  dumpTwoPathTable(g);
  for (i = 0; i < g->num_A_nodes; i++) {
    for (j = 0; j < g->num_A_nodes; j++) {
      if (i != j) {
	assert(GET_A2PATH_ENTRY(g, i, j) == twoPaths(g, i, j));
      }
    }
  }
  assert(g->num_A_nodes + g->num_B_nodes == g->num_nodes);
  for (i = g->num_A_nodes; i < g->num_A_nodes + g->num_B_nodes; i++) {
    for (j = g->num_A_nodes; j < g->num_A_nodes + g->num_B_nodes; j++) {
      if (i != j) {
	assert(GET_B2PATH_ENTRY(g, i, j) == twoPaths(g, i, j));
      }
    }
  }
#endif /*TWOPATH_LOOKUP*/

  for (i = 0; i < g->num_A_nodes; i++) {
    for (j = 0; j < g->num_B_nodes; j++) {
      if (i != j) {
        assert(isEdge(g, i, j) == isEdge(g, j, i));
      }
    }
  }

  /* just change stats (no changes to graph) */
  printf("testing change stats\n");
  gettimeofday(&start_timeval, NULL);
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
      fprintf(stderr, "edge %u -- %u already exists\n", i, j);
      continue;
    }

    printf("i = %d, j = %d, changeC4 = %g, changeKsp = %g, changeKsa = %g, changeKca = %g, changeKcp = %g, changeSa2 = %g, changeSp2 = %g, changeSa3 = %g, changeSp3 = %g, changeL3 = %g\n", i, j - g->num_A_nodes,
	   changeFourCycles(g, i, j, lambda), //XXX changeBipartiteFourCycle(g, i, j, lambda),
           changeBipartiteAltStarsA(g, i, j, lambda),
           changeBipartiteAltStarsB(g, i, j, lambda),
           changeBipartiteAltKCyclesB(g, i, j, lambda),
           changeBipartiteAltKCyclesA(g, i, j, lambda),
	   changeBipartiteTwoStarsB(g, i, j, lambda),
	   changeBipartiteTwoStarsA(g, i, j, lambda),
	   changeBipartiteThreeStarsB(g, i, j, lambda),
	   changeBipartiteThreeStarsA(g, i, j, lambda),
	   changeBipartiteThreePath(g, i, j, lambda));
    num_tests++;
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
  gettimeofday(&end_timeval, NULL);
  timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
  etime = 1000 * elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
  fprintf(stderr, "Change stats computations took %.2f s\n", (double)etime/1000);
  

  /* add edges and update graph and 2-path hash tables */
  printf("testing add edges\n");
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
      if (i >= g->num_nodes) {
        fprintf(stderr, "bad i node num %d\n", i);
        exit(1);
      }
      if (j >= g->num_nodes) {
        fprintf(stderr, "bad j node num %d\n", j);
        exit(1);
      }
      j += g->num_A_nodes; /* in file is 0 .. N-Np-1 but we number Np .. N-1 */      
    } else {
      do {
	i = rand() % g->num_A_nodes;
	j = g->num_A_nodes + rand() % g->num_B_nodes;
      } while (isEdge(g, i, j));
    }
    if (isEdge(g, i, j)) {
      continue;
    }
    insertEdge(g, i, j);
    /* insertEdge() called updateTwoPathsMatrices() itself */
#ifdef TWOPATH_LOOKUP
    printf("i = %d, j = %d, num_edges = %d, ", i, j - g->num_A_nodes, g->num_edges);
    dumpTwoPathTable(g);
#else
    printf("i = %d, j = %d, num_edges = %d\n", i, j - g->num_A_nodes, g->num_edges);
#endif /*TWOPATH_LOOKUP*/
    num_tests++;
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

  
  /* delete edges and update graph and 2-path hash tables */
  printf("testing delete edges\n");
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
      if (i >= g->num_nodes) {
        fprintf(stderr, "bad i node num %d\n", i);
        exit(1);
      }
      if (j >= g->num_nodes) {
        fprintf(stderr, "bad j node num %d\n", j);
        exit(1);
      }
      j += g->num_A_nodes; /* in file is 0 .. N-Np-1 but we number Np .. N-1 */      
    } else {
      do {
	i = rand() % g->num_A_nodes;
	j = g->num_A_nodes + rand() % g->num_B_nodes;
      } while (!isEdge(g, i, j));
    }
    if (!isEdge(g, i, j)) {
      continue;
    }
    removeEdge(g, i, j);
    /* removeEdge() calles updateTwoPathsMatrices() itself */
#ifdef TWOPATH_LOOKUP
    printf("i = %d, j = %d, num_edges = %d, ", i, j - g->num_A_nodes, g->num_edges);
    dumpTwoPathTable(g);
    #else
    printf("i = %d, j = %d, num_edges = %d\n", i, j - g->num_A_nodes, g->num_edges);
#endif /*TWOPATH_LOOKUP*/
    num_tests++;
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
