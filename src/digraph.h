#ifndef DIGRAPH_H
#define DIGRAPH_H
/*****************************************************************************
 * 
 * File:    digraph.h
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * Directed graph data structure. Stored as arc lists (both forward and
 * a "reversed" version, for fast iteration over both in- and out- neighbours)
 * and fast lookup hash tables for two-paths, and flat arcs list for fast 
 * finding of random arc.
 *
 * Nodes are numbered 0 .. n-1.
 *
 *
 ****************************************************************************/

#include <stdio.h>
#include "utils.h"
#include "uthash.h"

#define BIN_NA  -1  /* value for binary missing data (otherwise 0 or 1) */
#define CAT_NA  -1  /* value for catagorical missing data (otherwise >= 0) */

typedef struct nodepair_s /* pair of nodes (i, j) */
{
  uint_t  i;    /* from node */
  uint_t  j;    /* to node */
} nodepair_t;



/* combine 32 bit indices i and j into single 64 bit key for hash table */
#define MAKE_KEY64(i, j) (((uint64_t)(i) << 32) | ((j) & 0xffffffff))

/* uthash hash table entry hsa 64 bit key (32 bit i and j indices converted
   to 64 bit with MAKE_KEY64 macro) and 32 bit value (number of two-paths) */
typedef struct {
  uint64_t       key;   /* i, j indices packed into 64 bits withe MAKE_KEY64 */
  uint32_t       value; /* count of two-paths between i and j in key */
  UT_hash_handle hh;    /* uthash hash handle */
} twopath_record_t;

typedef struct digraph_s
{
  uint_t   num_nodes;  /* number of nodes */
  uint_t   num_arcs;   /* number of arcs */
  uint_t  *outdegree;  /* for each node, number of nodes it has an arc to */
  uint_t **arclist;    /* arc adjacency lists: for each node i, array of
                          outdegree[i] nodes it has an arc to */
  uint_t  *indegree;   /* for each node, number of nodes that have an arc to it*/
  uint_t **revarclist; /* reverse arc adjacency list: for each node i, array of 
                          indegree[i] nodes that have an arc to it */
  nodepair_t *allarcs; /* list of all arcs specified as i->j for each */

  /* the keys for hash tables are 64 bits: 32 bits each for i and j index */
  twopath_record_t *mixTwoPathHashTab; /* hash table counting two-paths */
  twopath_record_t *inTwoPathHashTab;  /* hash table counting in-two-paths */
  twopath_record_t *outTwoPathHashTab; /* hash table counting out-two-paths */

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

  /* use for GeoDistance, need to mark continuous attributes for lat/long */
  uint_t latitude_index;  /* index in digraph contattr of latitude */
  uint_t longitude_index; /* index in digraph contattr of longitude */

  /* snowball sampling information, only used for conditional estimation */
  uint_t *zone;        /* for each node, snowball sampling zone (0 for seeds) */
  uint_t max_zone;     /* highest zone number (zone number of outermost wave) */
  uint_t num_inner_nodes;/*number of nodes in inner waves (all but last zone)*/
  uint_t *inner_nodes; /* id of each of the num_inner_nodes inner wave nodes */
  uint_t *prev_wave_degree; /* for each  node, number of edges 
                               to/from a node in earlier wave (node zone -1 ) */
} digraph_t;

uint_t get_twopath_entry(twopath_record_t *h, uint_t i, uint_t j);
  
digraph_t *load_digraph_from_arclist_file(FILE *pajek_file,
                                          const char *binattr_filename,
                                          const char *catattr_filename,
                                          const char *contattr_filename);

double density(const digraph_t *g); /* graph density of g */
bool isArc(const digraph_t *g, uint_t i, uint_t j); /* test if arc i->j is in g */
bool isArcIgnoreDirection(const digraph_t *g, uint_t i, uint_t j); /* test if arc i->j or j->i is in g */

/* these two version do not update the allarcs flat arclist */
void insertArc(digraph_t *g, uint_t i, uint_t j); /* add arc i->j to g */
void removeArc(digraph_t *g, uint_t i, uint_t j); /* delete arc i->j from g */

/* this two versions update the allarcs flat arclist also */
void insertArc_allarcs(digraph_t *g, uint_t i, uint_t j); /* add arc i->j to g */
void removeArc_allarcs(digraph_t *g, uint_t i, uint_t j, uint_t arcidx); /* delete arc i->j from g */

digraph_t *allocate_digraph(uint_t num_vertices);
void free_digraph(digraph_t *g);
void dump_digraph_arclist(const digraph_t *g);
void print_data_summary(const digraph_t *g);
void print_zone_summary(const digraph_t *g);

void write_digraph_arclist_to_file(FILE *fp, const digraph_t *g);

int add_snowball_zones_to_digraph(digraph_t *g, const char *zone_filename);
void dump_zone_info(const digraph_t *g);

#endif /* DIGRAPH_H */

