#ifndef GRAPH_H
#define GRAPH_H
/*****************************************************************************
 * 
 * File:    graph.h
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * Directed or undirected graph data structure. Also handles bipartite
 * (two-mode)
 *
 * For directed, stored as arc lists (both forward and a "reversed"
 * version, for fast iteration over both in- and out- neighbours).
 * Also, fast lookup hash tables for two-paths, and flat arcs or edges
 * list for fast finding of random arc.
 *
 * Nodes are numbered 0 .. n-1.
 *
 *
 * Preprocessor defines used:
 *
 *    TWOPATH_LOOKUP      - use two-path lookup tables (arrays by default)
 *    TWOPATH_HASHTABLES  - use hash tables (only if TWOPATH_LOOKUP defined)
 *
 ****************************************************************************/

#include <stdio.h>
#include "utils.h"
#ifdef TWOPATH_LOOKUP
#ifdef TWOPATH_HASHTABLES
#include "uthash.h"
#endif /* TWOPATH_HASHTABLES */
#endif /* TWOPATH_LOOKUP */


#define BIN_NA  -1  /* value for binary missing data (otherwise 0 or 1) */
#define CAT_NA  -1  /* value for catagorical missing data (otherwise >= 0) */



/* set element type, each element in array is either present, absent, or NA */
typedef enum set_elem_e {
  SET_ELEM_NA        = -1,
  SET_ELEM_ABSENT    =  0,
  SET_ELEM_PRESENT   =  1
} set_elem_e;

typedef struct nodepair_s /* pair of nodes (i, j) */
{
  uint_t  i;    /* from node */
  uint_t  j;    /* to node */
} nodepair_t;


/* Mode (node type) of a node for bipartite (two-mode) networks */
typedef enum bipartite_node_mode_e {
  MODE_INVALID = -1,
  MODE_A      =  0,
  MODE_B       = 1
} bipartite_node_mode_e;

#ifdef TWOPATH_LOOKUP
#ifdef TWOPATH_HASHTABLES
/* uthash hash table entry has (i,j) as key and number of tw-paths as value */
typedef struct {
  nodepair_t     key;   /* i, j indices */
  uint32_t       value; /* count of two-paths between i and j in key */
  UT_hash_handle hh;    /* uthash hash handle */
} twopath_record_t;

#define GET_MIX2PATH_ENTRY(g, i, j) get_twopath_entry((g)->mixTwoPathHashTab, (i), (j))
#define GET_IN2PATH_ENTRY(g, i, j) get_twopath_entry((g)->inTwoPathHashTab, (i), (j))
#define GET_OUT2PATH_ENTRY(g, i, j) get_twopath_entry((g)->outTwoPathHashTab, (i), (j))
#define GET_2PATH_ENTRY(g, i, j) get_twopath_entry((g)->twoPathHashTab, (i), (j)) /* undirected */
#define GET_A2PATH_ENTRY(g, i, j) get_twopath_entry((g)->twoPathHashTabA, (i), (j)) /* bipartite */
#define GET_B2PATH_ENTRY(g, i, j) get_twopath_entry((g)->twoPathHashTabB, (i), (j)) /* bipartite */
#else
#define GET_MIX2PATH_ENTRY(g, i, j) ((g)->mixTwoPathMatrix[INDEX2D((i), (j), (g)->num_nodes)])
#define GET_IN2PATH_ENTRY(g, i, j) ((g)->inTwoPathMatrix[INDEX2D((i), (j), (g)->num_nodes)])
#define GET_OUT2PATH_ENTRY(g, i, j) ((g)->outTwoPathMatrix[INDEX2D((i), (j), (g)->num_nodes)])
#define GET_2PATH_ENTRY(g, i, j) ((g)->twoPathMatrix[INDEX2D((i), (j), (g)->num_nodes)]) /* undirected */
#define GET_A2PATH_ENTRY(g, i, j) ((g)->twoPathMatrixA[INDEX2D((i), (j), (g)->num_B_nodes)]) /* bipartite */
#define GET_B2PATH_ENTRY(g, i, j) ((g)->twoPathMatrixB[INDEX2D((i)-(g)->num_A_nodes, (j)-(g)->num_A_nodes, (g)->num_B_nodes)]) /* bipartite */ /* Note subtracting num_A_nodes as B nodes are numbered num_A_nodes .. num_nodes */
#endif /* TWOPATH_HASHTABLES */
#else /* not using two-path lookup tables (either arrays or hashtables) */
#define GET_MIX2PATH_ENTRY(g, i, j) mixTwoPaths((g), (i), (j))
#define GET_OUT2PATH_ENTRY(g, i, j) outTwoPaths((g), (i), (j))
#define GET_IN2PATH_ENTRY(g, i, j) inTwoPaths((g), (i), (j))
#define GET_2PATH_ENTRY(g, i, j) twoPaths((g), (i), (j)) /* undirected */
/* Note bipartite undirected can just use same twoPaths() function as one-mode*/
#define GET_A2PATH_ENTRY(g, i, j) twoPaths((g), (i), (j)) /* bipartite */
#define GET_B2PATH_ENTRY(g, i, j) twoPaths((g), (i), (j)) /* bipartite */
#endif /* TWOPATH_LOOKUP */

typedef struct graph_s
{
  uint_t   num_nodes;  /* number of nodes */

  bool     is_directed;  /* TRUE for directed graph, else undirected */
  bool     is_bipartite; /* TRUE for two-mode graph, else one-mode */

  /* 
   * Directed graph, used only if is_directed 
   */
  uint_t   num_arcs;   /* number of arcs */
  uint_t  *outdegree;  /* for each node, number of nodes it has an arc to */
  uint_t **arclist;    /* arc adjacency lists: for each node i, array of
                          outdegree[i] nodes it has an arc to */
  uint_t  *indegree;   /* for each node, number of nodes that have an arc to it*/
  uint_t **revarclist; /* reverse arc adjacency list: for each node i, array of 
                          indegree[i] nodes that have an arc to it */
  nodepair_t *allarcs; /* list of all arcs specified as i->j for each. */

#ifdef TWOPATH_LOOKUP
#ifdef TWOPATH_HASHTABLES
  /* the keys for hash tables are 64 bits: 32 bits each for i and j index */
  twopath_record_t *mixTwoPathHashTab; /* hash table counting two-paths */
  twopath_record_t *inTwoPathHashTab;  /* hash table counting in-two-paths */
  twopath_record_t *outTwoPathHashTab; /* hash table counting out-two-paths */
#else /* using array not hashtables for two-path lookup */
  uint_t *mixTwoPathMatrix; /* n x n contiguous matrix counting two-paths */
  uint_t *inTwoPathMatrix;  /* n x n contiguous matrix counting in-two-paths */
  uint_t *outTwoPathMatrix; /* n x n contiguous matrix counting out-two-paths */
#endif /*TWOPATH_HASHTABLES*/
#endif /* TWOPATH_LOOKUP */
  

  /*
   * Undirected graph, used only if !is_directed 
   */
  uint_t   num_edges;  /* number of edges */
  uint_t  *degree;     /* for each node, number of nodes it has an edge to */
  uint_t **edgelist;   /* edge adjacency lists: for each node i, array of
                          degree[i] nodes it has an arc to.
			  Note that for each edge i--j there will be two
			  entries, one for i and one for j */
  nodepair_t *alledges;/* list of all edges specified as i--j for each. */

#ifdef TWOPATH_LOOKUP
#ifdef TWOPATH_HASHTABLES
  /* the keys for hash tables are 64 bits: 32 bits each for i and j index */
  twopath_record_t *twoPathHashTab; /* hash table counting two-paths */
#else /* using array not hashtables for two-path lookup */
  uint_t *twoPathMatrix; /* n x n contiguous matrix counting two-paths */
#endif /*TWOPATH_HASHTABLES*/
#endif /* TWOPATH_LOOKUP */


  /*
   * Bipartite (two-mode) graph, used only if is_bipartite
   */

  /* In the two-mode (bipartite) graph, we refer to the two modes (node types)
     as A and B. There are only edges between nodes of different modes
     i.e. between a type A node and a type B node; there are no edges
     between nodes both of type A or both of type B.
     As usual the nodes are indexed 0 ... N-1, but for bipartite, all the
     type A nodes are first, followed by all the type B nodes, i.e.
     nodes 0 .. N_A are the A nodes and N_A .. N-1 are the type B nodes
     (where N = N_A + N_B).
  */
  uint_t num_A_nodes; /* number of mode A nodes. */
  uint_t num_B_nodes; /* number of mode B nodes.
			 num_B_nodes = num_nodes - num_A_nodes */

#ifdef TWOPATH_LOOKUP
#ifdef TWOPATH_HASHTABLES
  /* the keys for hash tables are 64 bits: 32 bits each for i and j index */
  twopath_record_t *twoPathHashTabA; /* hash table counting A-two-paths */
  twopath_record_t *twoPathHashTabB; /* hash table counting B-two-paths */
#else /* using array not hashtables for two-path lookup */
  uint_t *twoPathMatrixA; /* n_A x n_A contiguous matrix counting A-two-paths */
  uint_t *twoPathMatrixB; /* n_B x n_B contiguous matrix counting B-two-paths */
#endif /*TWOPATH_HASHTABLES*/
#endif /* TWOPATH_LOOKUP */


  /*
   * Fields used for both directed and undirected graphs
   */
  
  /* node attributes */
  uint_t   num_binattr;   /* number of binary attributes */
  char   **binattr_names; /* binary attribute names */
  int    **binattr;       /* binary attributes. For each binary attribute u,
                             binattr[u][i] is the value for node i or BIN_NA 
                             for missing data */
  uint_t   num_catattr;   /* number of categorical attributes */
  char   **catattr_names; /* categorical attributes names */
  int    **catattr;       /* categorical attribute. For each categorical
                             attribute u, catattr[u][i] is value for node i  or
                             CAT_NA for missing data */
  uint_t   num_contattr;  /* number of continuous attributes */
  char   **contattr_names;/* continuous attributes names */
  double **contattr;      /* continuous attribute. For each continuous
                             attribute u, contattr[u][i] is value for node i 
                             or IEEE NaN for missin data (test with isnan()) */
  uint_t        num_setattr;   /* number of set (of categorical) attributes */
  char       **setattr_names;  /* set attributes names */
  uint_t      *setattr_lengths;/* size of each set array cattr[u][i][] */
  set_elem_e ***setattr;       /* set attribute. For each attribute
                                  attribute u, catattr[u][i][j] is the
                                  present, absent, or NA value for
                                  element j of node i. */

  /* use for GeoDistance, need to mark continuous attributes for lat/long */
  uint_t latitude_index;  /* index in digraph contattr of latitude */
  uint_t longitude_index; /* index in digraph contattr of longitude */

  /* use for EuclidenDistance, need to mark contiuous attributes for x/y/z */
  uint_t x_index;         /* index in digraph contattr of x coordinate */
  uint_t y_index;         /* index in digraph contattr of y coordinate */
  uint_t z_index;         /* index in digraph contattr of z coordinate */

  /* 
   * snowball sampling information, only used for conditional estimation 
   */
  uint_t *zone;        /* for each node, snowball sampling zone (0 for seeds) */
  uint_t max_zone;     /* highest zone number (zone number of outermost wave) */
  uint_t num_inner_nodes;/*number of nodes in inner waves (all but last zone)*/
  uint_t *inner_nodes; /* id of each of the num_inner_nodes inner wave nodes */
  uint_t *prev_wave_degree; /* for each  node, number of edges 
                               to/from a node in earlier wave (node zone -1 ) */

  /* for digraphs only */
  uint_t num_inner_arcs;  /* number of arcs in inner waves, length of
                             allinnerarcs list */
  nodepair_t *allinnerarcs; /* list of all inner wave arcs specified
                             * as i->j for each. */
  /* for undirected graphs only */
  uint_t num_inner_edges;  /* number of edges in inner waves, length of
                             allinneredges list */
  nodepair_t *allinneredges; /* list of all inner wave edges specified
                             * as i--j for each. */

  /* 
   * term (time period) information, only used for citation ERGM (cERGM) 
   * which is only for directed graphs
   */
  uint_t *term;    /* for each node, the sequential time period (term) */
  uint_t max_term;  /* highest term number (time period) 0 ... max_term */
  uint_t num_maxterm_nodes; /* number of nodes in last (latest) term */
  uint_t *maxterm_nodes; /* id of each of the num_maxterm_nodes in the
			    latest term i.e. with term == max_term */
  uint_t num_maxtermsender_arcs;  /* number of arcs with sender in max. term,
				     length of all_maxtermsender_arcs list */
  nodepair_t *all_maxtermsender_arcs; /* list of all arcs where i->j where
					 i has maximum term value max_term,
					 specified as i->j for each. */
  
  
} graph_t;

#ifdef TWOPATH_LOOKUP
#ifdef TWOPATH_HASHTABLES
uint_t get_twopath_entry(twopath_record_t *h, uint_t i, uint_t j);
#endif /*TWOPATH_HASHTABLES */
#else /* not using two-path lookup tables (either arrays or hashtables) */
uint_t mixTwoPaths(const graph_t *g, uint_t i, uint_t j);
uint_t outTwoPaths(const graph_t *g, uint_t i, uint_t j);
uint_t inTwoPaths(const graph_t *g, uint_t i, uint_t j);
uint_t twoPaths(const graph_t *g, uint_t i, uint_t j);
#endif /*TWOPATH_LOOKUP */
  
double num_graph_dyads(const graph_t *g, bool allowLoops); /*max possible edges in g*/
double num_graph_inner_dyads(const graph_t *g); /*max edges beteen inner nodes*/
uint_t num_arcs_or_edges(const graph_t *g); /* number of arcs or edges in g */
uint_t num_inner_arcs_or_edges(const graph_t *g);
double density(const graph_t *g, bool allowLoops); /* graph density of g */
bool isArc(const graph_t *g, uint_t i, uint_t j); /* test if arc i->j is in g */
bool isEdge(const graph_t *g, uint_t i, uint_t j); /* test if edge i--j in in g */
bool isArcIgnoreDirection(const graph_t *g, uint_t i, uint_t j); /* test if arc i->j or j->i is in g */
bool isArcOrEdge(const graph_t *g, uint_t i, uint_t j); /* test if arc i->j is in g if directed or edge i -- j if undirected */

/* these two version do not update the allarcs flat arclist */
void insertArc(graph_t *g, uint_t i, uint_t j); /* add arc i->j to g */
void removeArc(graph_t *g, uint_t i, uint_t j); /* delete arc i->j from g */

/* these two version do not update the alledges flat edgelist */
void insertEdge(graph_t *g, uint_t i, uint_t j); /* add edge i--j to g */
void removeEdge(graph_t *g, uint_t i, uint_t j); /* delete edge i--j from g */

/* this two versions update the allarcs flat arclist also */
void insertArc_allarcs(graph_t *g, uint_t i, uint_t j); /* add arc i->j to g */
void removeArc_allarcs(graph_t *g, uint_t i, uint_t j, uint_t arcidx); /* delete arc i->j from g */

/* this two versions update the alledges flat edgelist also */
void insertEdge_alledges(graph_t *g, uint_t i, uint_t j); /* add edge i--j to g */
void removeEdge_alledges(graph_t *g, uint_t i, uint_t j, uint_t edgeidx); /* delete edge i--j from g */

/* this two versions update the allinnerarcs flat arclist also */
void insertArc_allinnerarcs(graph_t *g, uint_t i, uint_t j); /* add arc i->j to g */
void removeArc_allinnerarcs(graph_t *g, uint_t i, uint_t j, uint_t arcidx); /* delete arc i->j from g */

/* this two versions update the allinneredges flat edgelist also */
void insertEdge_allinneredges(graph_t *g, uint_t i, uint_t j); /* add edge i--j to g */
void removeEdge_allinneredges(graph_t *g, uint_t i, uint_t j, uint_t edgeidx); /* delete edge i--j from g */

/* this two versions update the all_maxtermsender_arcs flat arclist also */
void insertArc_all_maxtermsender_arcs(graph_t *g, uint_t i, uint_t j); /* add arc i->j to g */
void removeArc_all_maxtermsender_arcs(graph_t *g, uint_t i, uint_t j, uint_t arcidx); /* delete arc i->j from g */


graph_t *allocate_graph(uint_t num_vertices, bool is_directed, bool is_bipartite);
void free_graph(graph_t *g);
void dump_graph_arclist(const graph_t *g);
void print_data_summary(const graph_t *g, bool allowLoops);
void print_zone_summary(const graph_t *g);

void write_graph_arclist_to_file(FILE *fp, const graph_t *g);

int add_snowball_zones_to_graph(graph_t *g, const char *zone_filename);
void dump_zone_info(const graph_t *g);

int parse_category_set(char *str, bool firstpass, uint_t *size,
                       set_elem_e *setval);

int load_attributes(graph_t *g, 
                    const char *binattr_filename,
                    const char *catattr_filename,
                    const char *contattr_filename,
                    const char *setattr_filename);

uint_t get_num_vertices_from_arclist_file(FILE *pajek_file);

int add_cergm_terms_to_digraph(graph_t *g, const char *term_filename);
void dump_term_info(const graph_t *g);
void print_term_summary(const graph_t *g);
uint_t num_loops(const graph_t *g);
bool has_loop(const graph_t *t, uint_t u);

/* these versions insert/remove arc or edge as appropriate for g */
void insertArcOrEdge(graph_t *g, uint_t i, uint_t j);
void removeArcOrEdge(graph_t *g, uint_t i, uint_t j);
void insertArcOrEdge_updatelist(graph_t *g, uint_t i, uint_t j);
void removeArcOrEdge_updatelist(graph_t *g, uint_t i, uint_t j, uint_t idx);
void insertArcOrEdge_updateinnerlist(graph_t *g, uint_t i, uint_t j);
void removeArcOrEdge_updateinnerlist(graph_t *g, uint_t i, uint_t j, uint_t idx);

bipartite_node_mode_e bipartite_node_mode(const graph_t *g, uint_t i);

#endif /* GRAPH_H */

