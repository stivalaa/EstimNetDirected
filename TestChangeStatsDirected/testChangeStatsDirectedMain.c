/*****************************************************************************
 * 
 * File:    testChangeStatsDirected.c
 * Author:  Alex Stivala
 * Created: October 2017
 *
 * Test directed change stats and two-path table update.
 *
 *
 * Usage:  testChangeStatsDirected  <in_edgelistfile> [nodenums]
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
#include "changeStatisticsDirected.h"
#include "loadDigraph.h"

#define DEFAULT_NUM_TESTS 1000

#ifdef TWOPATH_LOOKUP
/* get stats and dump mix-two-path hash table */
static void dumpTwoPathTables(const graph_t *g) {
  uint_t inSum,outSum,mixSum,inMax,outMax,mixMax;
  uint_t inNnz, outNnz, mixNnz;
#ifdef TWOPATH_HASHTABLES
  twopath_record_t *r;
  uint_t val;
#else
  uint_t i, j;
#endif /*TWOPATH_HASHTABLES*/
  
  inSum = outSum = mixSum = 0;
  inMax = outMax = mixMax = 0;
  inNnz = outNnz = mixNnz = 0;

#ifdef TWOPATH_HASHTABLES
  mixNnz = HASH_COUNT(g->mixTwoPathHashTab);
  inNnz = HASH_COUNT(g->inTwoPathHashTab);
  outNnz = HASH_COUNT(g->outTwoPathHashTab);

  for (r = g->mixTwoPathHashTab; r !=NULL; r = r->hh.next) {
    val = r->value;
    mixSum += val;
    if (val > mixMax) {
      mixMax = val;
    }
  }
  for (r = g->inTwoPathHashTab; r !=NULL; r = r->hh.next) {
    val = r->value;
    inSum += val;
    if (val > inMax) {
      inMax = val;
    }
  }
  for (r = g->outTwoPathHashTab; r !=NULL; r = r->hh.next) {
    val = r->value;
    outSum += val;
    if (val > outMax) {
      outMax = val;
    }
  }
#else
 for (i = 0; i < g->num_nodes; i++) {
    for (j = 0; j < g->num_nodes; j++) {
      mixSum += g->mixTwoPathMatrix[INDEX2D(i, j, g->num_nodes)];
      inSum += g->inTwoPathMatrix[INDEX2D(i, j, g->num_nodes)];
      outSum += g->outTwoPathMatrix[INDEX2D(i, j, g->num_nodes)];
      if (g->mixTwoPathMatrix[INDEX2D(i, j, g->num_nodes)] > 0) {
        mixNnz++;
      }
      if (g->mixTwoPathMatrix[INDEX2D(i, j, g->num_nodes)] > mixMax) {
        mixMax = g->mixTwoPathMatrix[INDEX2D(i, j, g->num_nodes)];
      }
      if (g->inTwoPathMatrix[INDEX2D(i, j, g->num_nodes)] > 0) {
        inNnz++;
      }
      if (g->inTwoPathMatrix[INDEX2D(i, j, g->num_nodes)] > inMax) {
        inMax = g->inTwoPathMatrix[INDEX2D(i, j, g->num_nodes)];
      }
      if (g->outTwoPathMatrix[INDEX2D(i, j, g->num_nodes)] > 0) {
        outNnz++;
      }
      if (g->outTwoPathMatrix[INDEX2D(i, j, g->num_nodes)] > outMax) {
        outMax = g->outTwoPathMatrix[INDEX2D(i, j, g->num_nodes)];
      }
    }
  }
 #endif /*TWOPATH_HASHTABLES*/

  printf("mix2p sum = %u, max = %u\n", mixSum, mixMax);
  printf("in2p sum = %u, max = %u\n", inSum, inMax);
  printf("out2p sum = %u, max = %u\n", outSum, outMax);
  printf("mix nnz = %u (%.4f%%)\n", mixNnz,
         100*(double)mixNnz/(g->num_nodes*g->num_nodes));
  printf("in nnz = %u (%.4f%%)\n", inNnz,
         100*(double)inNnz/(g->num_nodes*g->num_nodes));
  printf("out nnz = %u (%.4f%%)\n", outNnz,
         100*(double)outNnz/(g->num_nodes*g->num_nodes));
}
#endif /*TWOPATH_LOOKUP*/


int main(int argc, char *argv[]) 
{
  const double lambda = 2.0;  /* always use decay lambda = 2.0 on tests */
  char buf[1024];
  uint_t i,j;
  char *arclist_filename = NULL;
  FILE *file           = NULL;
  uint_t     num_nodes = 0;
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
  arclist_filename = argv[1];
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

  if (!(file = fopen(arclist_filename, "r"))) {
    fprintf(stderr, "error opening file %s (%s)\n", 
            arclist_filename, strerror(errno));
    return -1;
  }
  
  gettimeofday(&start_timeval, NULL);
#ifdef TWOPATH_LOOKUP
  fprintf(stderr, "loading arc list from %s and building two-path tables...",
         arclist_filename);
#else
    fprintf(stderr, "loading arc list from %s...",
         arclist_filename);
#endif
  num_nodes = get_num_vertices_from_arclist_file(file); /* closes file */
  g = allocate_graph(num_nodes, TRUE);
  if (!(file = fopen(arclist_filename, "r"))) {
    fprintf(stderr, "error opening file %s (%s)\n", 
            arclist_filename, strerror(errno));
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
  dumpTwoPathTables(g);
#endif /*TWOPATH_LOOKUP*/

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
      if (i >= g->num_nodes) {
        fprintf(stderr, "bad i node num %d\n", i);
        exit(1);
      }
      if (j >= g->num_nodes) {
        fprintf(stderr, "bad j node num %d\n", j);
        exit(1);
      }
    } else {
      i = rand() % g->num_nodes;
      j = rand() % g->num_nodes;
    }
    if (i == j) {
      continue;
    }
    printf("i = %d, j = %d, changeOutKStars = %g, changeInKStars = %g, changeDiTKTriangles = %g, changeA2pTD = %g, changeDiCKTriangles = %g, changeDiUKTriangles = %g, changeDiDKTriangles = %g, changeDiUAltTwoPaths = %g, changeSource = %g, changeSink = %g, changeDiIso = %g, changeTwoMixStar = %g, change030c = %g, change030t = %g, changeIn2star = %g, changeOut2star = %g\n", i, j,
           changeAltOutStars(g, i, j, lambda),
           changeAltInStars(g, i, j, lambda),
           changeAltKTrianglesT(g, i, j, lambda),
           changeAltTwoPathsTD(g, i, j, lambda),
           changeAltKTrianglesC(g, i, j, lambda),
           changeAltKTrianglesU(g, i, j, lambda),
           changeAltKTrianglesD(g, i, j, lambda),
           changeAltTwoPathsU(g, i, j, lambda),
           changeSource(g, i, j, lambda),
           changeSink(g, i, j, lambda),
           changeIsolates(g, i, j, lambda),
	   changeTwoPath(g, i, j, lambda),
	   changeCyclicTriad(g, i, j, lambda),
	   changeTransitiveTriad(g, i, j, lambda),
	   changeInTwoStars(g, i, j, lambda),
	   changeOutTwoStars(g, i, j, lambda)
      );
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
  

  /* add arcs and update graph and 2-path hash tables */
  printf("testing add arcs\n");
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
    } else {
      do {
        i = rand() % g->num_nodes;
        j = rand() % g->num_nodes;
      } while (i == j || isArc(g, i, j));
    }
    if (i == j || isArc(g, i, j)) {
      continue;
    }
    insertArc(g, i, j);
    /* insertArc() called updateTwoPathsMatrices() itself */
#ifdef TWOPATH_LOOKUP
    printf("i = %d, j = %d, num_arcs = %d, ", i, j, g->num_arcs);
    dumpTwoPathTables(g);
#else
    printf("i = %d, j = %d, num_arcs = %d\n", i, j, g->num_arcs);
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

  
  /* delete arcs and update graph and 2-path hash tables */
  printf("testing delete arcs\n");
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
    } else {
      do {
        i = rand() % g->num_nodes;
        j = rand() % g->num_nodes;
      } while (i == j || !isArc(g, i, j));
    }
    if (i == j || !isArc(g, i, j)) {
      continue;
    }
    removeArc(g, i, j);
    /* removeArc() calles updateTwoPathsMatrices() itself */
#ifdef TWOPATH_LOOKUP
    printf("i = %d, j = %d, num_arcs = %d, ", i, j, g->num_arcs);
    dumpTwoPathTables(g);
    #else
    printf("i = %d, j = %d, num_arcs = %d\n", i, j, g->num_arcs);
#endif /*TWOPATH_LOOKUP*/
    num_tests++;
    if (!readNodeNums && num_tests >= DEFAULT_NUM_TESTS) {
      break;
    }
  }
  
  if (readNodeNums) {
    fclose(nodenumfile);
  }
  exit(0);
}
