/*****************************************************************************
 * 
 * File:    loadGraph.c
 * Author:  Alex Stivala
 * Created: October 2017
 *
 * Load digraph or graph from Pajek format arc list file and optionally compute
 * statistics corresponding to ERGM parameters.
 *
 * Preprocessor defines used:
 *
 *    TWOPATH_LOOKUP      - use two-path lookup tables (arrays by default)
 *    TWOPATH_HASHTABLES  - use hash tables (only if TWOPATH_LOOKUP defined)
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <ctype.h>
#include "utils.h"
#include "loadGraph.h"

   
/*****************************************************************************
 *
 * local constants
 *
 ****************************************************************************/

static const size_t BUFSIZE = 16384;  /* line buffer size for reading files */

/*****************************************************************************
 *
 * externally visible functions
 *
 ****************************************************************************/


/*
 * Build graph or digraph from Pajek format arc list file.
 * Note this builds the graph structure only, not the node attributes.
 * These are read and built separately by calling load_attributes().
 * The graph g must be already allocated by allocate_graph() with the
 * correct number of nodes, which can be found by
 * get_num_nodes(pajek_file).
 * The pre-allocated graph g also must have is_directed correctly
 * set for directed or undirected graph.
 *
 * In the Pajek format *vertices at top, then followed by one line for each
 * vertex (just vertex number) then *arcs (or *edges for undirected)
 * followed by arcs list one per
 * line. In this program the nodes must be numbered 1..N.
 *
 * For bipartite (two-mode) networks,
 * the first lines should be e.g.
 * *vertices 36 10
 * first number is total number of nodes
 * second number is number of mode A nodes
 * the rest are mode B - conventionally in the affiliation
 * matrix the rows are mode A and the columns mode B, e.g. mode A is
 * actors and mode B is their affiliations.
 * They must be numbered 1 ... N where N = num_A + num_B
 * so nodes 1 .. num_A are type A and num_A+1 .. N are type B
 * see e.g. http://www.pfeffer.at/txt2pajek/txt2pajek.pdf
 *
 * Parameters:
 *   pajek_file   - Pajek format arclist file handle (open read).
 *                   Closed by this function at end.
 *   g             - (in/out) digraph object already allocated as above.
 *   computeStats  - if TRUE the observed graph sufficient statistics
 *                   are computed by accumulating the change statistics as
 *                   each arc is added to the digraph. The following parameters
 *                   are then used (if FALSE then they can all be passed in as
 *                   0 or NULL):
 *   n      - number of parameters (length of theta vector and total
 *            number of change statistic functions)
 *   n_attr - number of attribute change stats functions
 *   n_dyadic -number of dyadic covariate change stats funcs
 *   n_attr_interaction - number of attribute interaction change stats funcs
 *   change_stats_funcs - array of pointers to change statistics functions
 *                        length is n-n_attr-n_dyadic-n_attr_interaction
 *   lambda_values      - array of lambda values for change stats funcs
 *                        same length as change_stats_funcs
 *   attr_change_stats_funcs - array of pointers to change statistics functions
 *                             length is n_attr
 *   dyadic_change_stats_funcs - array of pointers to dyadic change stats funcs
 *                             length is n_dyadic
 *   attr_interaction_change_stats_funcs - array of pointers to attribute
 *                           interaction (pair) change statistics functions.
 *                           length is n_attr_interaction.
 *   attr_indices   - array of n_attr attribute indices (index into g->binattr
 *                    or g->catattr) corresponding to attr_change_stats_funcs
 *                    E.g. for Sender effect on the first binary attribute,
 *                    attr_indices[x] = 0 and attr_change_stats_funcs[x] =
 *                    changeSender
 *   attr_interaction_pair_indices - array of n_attr_interaction pairs
 *                          of attribute inidices similar to above but
 *                          for attr_interaction_change_setats_funcs which
 *                          requires pairs of indices.
 *   addChangeStats           - (Out) vector of n change stats for add moves
 *                              Allocated by caller. Assuming we start
 *                              with the empty graph, this will then
 *                              have the counts of statistics in the 
 *                              graph at the end.
 *  theta                     - parameter values (required by 
 *                              calcChangeStats for total but value
 *                              not used here)
 *
 * Return value:
 *    digraph object built from files (same as parameter g)
 *
 * Note this function calls exit() on error.
 */
graph_t *load_graph_from_arclist_file(FILE *pajek_file, graph_t *g,
                                          bool computeStats,
                                          uint_t n, uint_t n_attr,
                                          uint_t n_dyadic,
                                          uint_t n_attr_interaction,
                                          change_stats_func_t
                                                     *change_stats_funcs[],
                                          double lambda_values[],
                                          attr_change_stats_func_t
                                          *attr_change_stats_funcs[],
                                          dyadic_change_stats_func_t
                                          *dyadic_change_stats_funcs[],
                                          attr_interaction_change_stats_func_t
                                         *attr_interaction_change_stats_funcs[],
                                          uint_t attr_indices[],
                                          uint_pair_t
                                          attr_interaction_pair_indices[],
                                          double *addChangeStats,
                                          double theta[])
{
  int i, j;
  char *p;
  char *saveptr   = NULL; /* for strtok_r() */
  char *token     = NULL; /* from strtok_r() */
  const char *delims = " \t\r\n"; /* strtok_r() delimiters for header lines */
  int num_vertices = 0;
  int num_A_vertices = 0;
  uint_t l;
  double *changestats = NULL;
#ifdef DEBUG_MEMUSAGE
  uint_t k, total_degree = 0;
#ifdef TWOPATH_HASHTABLES
  struct timeval start_timeval, end_timeval, elapsed_timeval;
  int            etime;
#endif /*TWOPATH_HASHABLES */
#endif /* DEBUG_MEMUSAGE */
  const char *directed_arcs_start_string = "*arcs"; /* case insensitive */
  const char *undirected_edges_start_string = "*edges";
  const char *edges_start_string = g->is_directed ? directed_arcs_start_string :
    undirected_edges_start_string;

  char buf[BUFSIZE];
  /* the first lines should be e.g.
   * *vertices 36
   * for one mode, or
   * *vertices 36 10
   * for two-mode
   * for Pajek format
   */
  fgets(buf, sizeof(buf)-1, pajek_file);
  for (p = buf; *p !='\0'; p++) {
    *p = tolower(*p);
  }
  if (g->is_bipartite) {
    if (sscanf(buf, "*vertices %d %d\n", &num_vertices, &num_A_vertices) != 2) {
      fprintf(stderr, "ERROR: expected *vertices n n_A line but didn't find it\n");
      exit(1);
    }
  } else {
    if (sscanf(buf, "*vertices %d\n", &num_vertices) != 1) {
      fprintf(stderr, "ERROR: expected *vertices n line but didn't find it\n");
      exit(1);
    }
  }

  if ((uint_t)num_vertices != g->num_nodes) {
    fprintf(stderr, "ERROR: expected %u vertices but found %d\n",
            g->num_nodes, num_vertices);
    exit(1);
  }
  if (g->is_bipartite && (uint_t)num_A_vertices != g->num_A_nodes) {
    fprintf(stderr, "ERROR: expected %u vertices in mode A for two-mode network but found %d\n",
            g->num_A_nodes, num_A_vertices);
    exit(1);
  }
  
  do {
    fgets(buf, sizeof(buf)-1, pajek_file);
    rstrip(buf);
  } while (!feof(pajek_file) && strcasecmp(buf, edges_start_string) != 0);
  if (feof(pajek_file)) {
    fprintf(stderr, "did not find %s line\n", edges_start_string);
    exit(1);
  }
  if (!fgets(buf, sizeof(buf)-1, pajek_file)) {
    fprintf(stderr, "ERROR: attempting to read first arc  (%s)\n", strerror(errno));
    fclose(pajek_file);
    exit(1);
  }

  if (computeStats)
    changestats = (double *)safe_malloc(n*sizeof(double));
  
#ifdef TWOPATH_HASHTABLES
#ifdef DEBUG_MEMUSAGE
  gettimeofday(&start_timeval, NULL);
#endif /* DEBUG_MEMUSAGE */
#endif /*TWOPATH_HASHTABLES*/
  while (!feof(pajek_file)) {
    i = j = 0;
    token = strtok_r(buf, delims, &saveptr);
    if (!token)
      break; /* end on blank line (ignore rest of file, used for stats etc.) */
    if (sscanf(token, "%d", &i) != 1) {
      fprintf(stderr, "ERROR: bad arc start node %s\n", token ? token : "(null)");
      exit(1);
    }
    token = strtok_r(NULL, delims, &saveptr);
    if (!token || sscanf(token, "%d", &j) != 1) {
      fprintf(stderr, "ERROR: bad arc end node %s\n", token ? token : "(null)");
      exit(1);
    }
    token = strtok_r(NULL, delims, &saveptr);
    if (token) {
      printf("(warning) ignoring Pajek arc weight %s on edge (%d,%d)\n", token, i, j);
    }

    
    if (i < 1 || j < 1) {
      fprintf(stderr, "ERROR: node numbers start at 1, got %d,%d\n", i, j);
      exit(1);
    }
    if (i > num_vertices || j > num_vertices) {
      fprintf(stderr, "ERROR num vertices %d but got edge %d,%d\n", num_vertices, i, j);
      exit(1); 
    }

    i--; j--; /* convert to 0-based */

    if (g->is_bipartite) {
      if (bipartite_node_mode(g, i) == bipartite_node_mode(g, j)) {
        fprintf(stderr, "ERROR: network is two-mode but edge %d,%d is between two nodes of the same mode\n", i, j);
      }
    }

    if (computeStats) {
      /* accumulate change statistics in addChangeStats array */
      (void)calcChangeStats(g, i, j, n, n_attr, n_dyadic, n_attr_interaction,
                            change_stats_funcs,
                            lambda_values,
                            attr_change_stats_funcs,
                            dyadic_change_stats_funcs,
                            attr_interaction_change_stats_funcs,
                            attr_indices,
                            attr_interaction_pair_indices,
                            theta, FALSE, changestats);
      for (l = 0; l < n; l++)
        addChangeStats[l] += changestats[l];
    }

    
    /* TODO calling isArc all the time is inefficient: since we are using
       hash table anyway would be better to check these in temporary hash
       table structure (or build graph as hash table then convert to the
       adjacency list structures). */
    if (g->is_directed) {
      if (!isArc(g, i, j)){
        insertArc_allarcs(g, i, j); /* also update flat arclist allarcs */
      }
    } else {
      /* undirected */
      if (!isEdge(g, i, j)){
        insertEdge_alledges(g, i, j); /* also update flat edgelist alledges */
      }
    }

#ifdef TWOPATH_HASHTABLES
#ifdef DEBUG_MEMUSAGE
    if (g->is_directed && g->num_arcs % 1000 == 0){
      gettimeofday(&end_timeval, NULL);
      timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
      etime = 1000 * elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
      MEMUSAGE_DEBUG_PRINT(("%u arcs, %u mixTwoPathHashtTab entries, %u inTwoPathHashTab entries, %u outTwoPathHashTab entries (%.2f s)...\n",
                            g->num_arcs,
                            HASH_COUNT(g->mixTwoPathHashTab),
                            HASH_COUNT(g->inTwoPathHashTab),
                            HASH_COUNT(g->outTwoPathHashTab),
                            (double)etime/1000));
    }
#endif /* DEBUG_MEMUSAGE */
#endif /* TWOPATH_HASHTABLES */
    
    saveptr = NULL; /* reset strtok() for next line */
    if (!fgets(buf, sizeof(buf)-1, pajek_file)) {
      if (!feof(pajek_file)) {
        fprintf(stderr, "ERROR: attempting to read edge (%s)\n", strerror(errno));
        fclose(pajek_file);
        exit(1);
      }
    }
  }
  fclose(pajek_file);

#ifdef DEBUG_MEMUSAGE
  if (g->is_directed) {
    for (k = 0; k < g->num_nodes; k++) {
      total_degree += g->outdegree[k];
    }
    MEMUSAGE_DEBUG_PRINT(("Allocated additional %f MB (twice) for %u arcs\n",
                          ((double)sizeof(uint_t) * total_degree) / (1024*1024),
                          g->num_arcs));
#ifdef TWOPATH_HASHTABLES
    MEMUSAGE_DEBUG_PRINT(("MixTwoPath hash table has %u entries (%f MB) which is %f%% nonzero in dense matrix (which would have taken %f MB)\n",
                          HASH_COUNT(g->mixTwoPathHashTab),
                          ((double)HASH_COUNT(g->mixTwoPathHashTab)*2*
                           sizeof(twopath_record_t))/(1024*1024),
                          100*(double)HASH_COUNT(g->mixTwoPathHashTab) /
                          ((double)g->num_nodes*g->num_nodes),
                          ((double)sizeof(uint_t)*num_vertices*num_vertices) /
                          (1024*1024)));
    MEMUSAGE_DEBUG_PRINT(("MixTwoPath hash table overhead %f MB\n", (double)HASH_OVERHEAD(hh, g->mixTwoPathHashTab)/(1024*1024)));
    MEMUSAGE_DEBUG_PRINT(("InTwoPath hash table has %u entries (%f MB) which is %f%% nonzero in dense matrix (which would have taken %f MB)\n",
                          HASH_COUNT(g->inTwoPathHashTab),
                          ((double)HASH_COUNT(g->inTwoPathHashTab)*
                           (sizeof(twopath_record_t)))/(1024*1024),
                          100*(double)HASH_COUNT(g->inTwoPathHashTab) /
                          ((double)g->num_nodes*g->num_nodes),
                          ((double)sizeof(uint_t)*num_vertices*num_vertices) /
                          (1024*1024)));
    MEMUSAGE_DEBUG_PRINT(("InTwoPath hash table overhead %f MB\n", (double)HASH_OVERHEAD(hh, g->inTwoPathHashTab)/(1024*1024)));
    MEMUSAGE_DEBUG_PRINT(("OutTwoPath hash table has %u entries (%f MB) which is %f%% nonzero in dense matrix (which would have taken %f MB)\n",
                          HASH_COUNT(g->outTwoPathHashTab),
                          ((double)HASH_COUNT(g->outTwoPathHashTab)*2*
                           sizeof(twopath_record_t))/(1024*1024),
                          100*(double)HASH_COUNT(g->outTwoPathHashTab) /
                          ((double)g->num_nodes*g->num_nodes),
                          ((double)sizeof(uint_t)*num_vertices*num_vertices) /
                          (1024*1024)));
    MEMUSAGE_DEBUG_PRINT(("OutTwoPath hash table overhead %f MB\n", (double)HASH_OVERHEAD(hh, g->outTwoPathHashTab)/(1024*1024)));
#endif /* TWOPATH_HASHTABLES */
  }
#endif /* DEBUG_MEMUSAGE */

  free(changestats);

  return(g);
}
