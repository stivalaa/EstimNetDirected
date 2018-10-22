/*****************************************************************************
 * 
 * File:    digraph.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * Directed graph data structure. Stored as arc lists (both forward and
 * a "reversed" version, for fast iteration over both in- and out- neighbours)
 * and fast lookup matrices for two-paths, and also flat arcs list for fast
 * selection of an arc uniformly at random.
 *
 *
 ****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "digraph.h"

   
/*****************************************************************************
 *
 * constants
 *
 ****************************************************************************/

static const size_t BUFSIZE = 16384;  /* line buffer size for reading files */

static const char *NA_STRING = "NA"; /* string in attributes files to indicate
                                        missing data (case insensitive) */

   
/*****************************************************************************
 *
 * local functions
 *
 ****************************************************************************/

/*
 * Update entry for (i, j) in hashtable.
 *
 * Used for the two-path hash tables to count two-paths between nodes
 * i and j. Much more efficient to use a hash table (or some other way
 * of storing a sparse matrix) as typically only on the order of 1% to
 * 10% of entries are nonzero. A hashtable is good as we want quick
 * lookup but not necessarily iteration (for which CSR etc. might be
 * better).
 *
 * Parameters:
 *     h - hash table
 *     i - node id of source
 *     j - node id of destination
 *     incval - value to add to existing value (or insert if not exists)
 *
 * Return value:
 *     None.
 */
static void update_twopath_entry(khash_t(m64) *h, uint_t i, uint_t j,
                                 uint_t incval)
{
  int      absent, is_missing;
  uint64_t kiter;
  uint64_t key = ((uint64_t)i << 32) | (j & 0xffffffff);
  kiter = kh_get(m64, h, key);
  is_missing = (kiter == kh_end(h));
  if (is_missing) {
    kiter = kh_put(m64, h, key, &absent);
    if (absent == -1) {
      fprintf(stderr, "ERROR: hash table key insert failed\n");
      exit(-1);
    }
    assert(absent > 0); /* key was not present, tested with kiter above */
    kh_value(h, kiter) = incval;
  } else {
    kh_value(h, kiter) = kh_value(h, kiter) + incval;
  }
} 


/*
 * Update the two-paths matrices used for fast computation of change
 * statistics for either adding or removing arc i->j
 * The matrices in the digraph are updated in-place
 *
 * Parameters:
 *   g     - digraph
 *   i     - node arc is from
 *   j     - node arc is to
 *   isAdd - TRUE for inserting arc, FALSE for deleting arc
 *
 * Return value:
 *   None.
 */
static void updateTwoPathsMatrices(digraph_t *g, uint_t i, uint_t j, bool isAdd)
{
  uint_t v,k;
  int incval = isAdd ? 1 : -1;
  /* TODO change dense matrices to sparse (hash table or CSR etc.) for scalabiity */
  for (k = 0; k < g->outdegree[i]; k++) {
    v = g->arclist[i][k];
    if (v == i || v == j)
      continue;
    assert(isArc(g,i,v));
    g->outTwoPathMatrix[INDEX2D(v, j, g->num_nodes)] += incval;
    g->outTwoPathMatrix[INDEX2D(j, v, g->num_nodes)] += incval;
  }
  for (k = 0; k < g->indegree[j]; k++) {
    v = g->revarclist[j][k];
    if (v == i || v == j)
      continue;
    assert(isArc(g,v,j));
    g->inTwoPathMatrix[INDEX2D(v, i, g->num_nodes)] += incval;
    g->inTwoPathMatrix[INDEX2D(i, v, g->num_nodes)] += incval;
  }
  for (k = 0; k < g->indegree[i]; k++)  {
    v = g->revarclist[i][k];
    if (v == i || v == j)
      continue;
    assert(isArc(g,v,i));
    g->mixTwoPathMatrix[INDEX2D(v, j, g->num_nodes)]+=incval;
    update_twopath_entry(g->mixTwoPathHashTab, v, j, incval);

  }
  for (k = 0; k < g->outdegree[j]; k++) {
    v = g->arclist[j][k];
    if (v == i || v == j)
      continue;
    assert(isArc(g,j,v));
    g->mixTwoPathMatrix[INDEX2D(i, v, g->num_nodes)] += incval;
    update_twopath_entry(g->mixTwoPathHashTab, i, v, incval);
  }
}


/*
 * Load integer (binary or categorical) attributes from file.
 * The format of the file is a header line with whitespace
 * delimited attribute names, and each subsequent line
 * the attribute values for each attribute.
 * The first line (after the header) has the values for
 * node 0, then the next line node 1, and so on.
 * 
 * E.g.:
 *
 * gender class
 * 0      1
 * 1      2
 * 1      3
 *
 *
 * Valid values are integer >= 0 for categorical or 0 or 1 for binary, 
 * or NA (case insensitve) for missing data.
 *
 * Parameters:
 *   attr_filenname - filename of file to read
 *   num_nodes - number of nodes (must be this many values)
 *   isBinary  - TRUE if binary (only 0 or 1 allowed)
 *               also any integer is allowed
 *   out_attr_names - (Out) attribute names array
 *   out_attr_values - (Out) (*attr_values)[u][i] is value of attr u for node i
 * 
 * Return value:
 *   Number of attributes, or -1 on error.
 *
 * The attribute names and values arrays are allocated by 
 * this function.
 */
static int load_integer_attributes(const char *attr_filename,
                                   uint_t num_nodes,
                                   bool isBinary,
                                   char ***out_attr_names,
                                   int  ***out_attr_values)
{
  const char *delims    = " \t\r\n"; /* strtok_r() delimiters  */
  uint_t nodenum        = 0;   /* node number values are for */
  uint_t num_attributes = 0;   /* number of different attributes */
  uint_t thisline_values= 0;   /* number values read this line */
  char  **attr_names   = NULL; /* array of attribute names */
  int   **attr_values = NULL; /* attr_values[u][i] is value of attr u for node i */
  char *saveptr        = NULL; /* for strtok_r() */
  char *token          = NULL; /* from strtok_r() */
  FILE *attr_file;
  char buf[BUFSIZE];
  uint_t  i;
  int     val;

  if (!(attr_file = fopen(attr_filename, "r"))) {
    fprintf(stderr, "ERROR: could not open attribute file %s (%s)\n",
            attr_filename, strerror(errno));
    return -1;
  }
  if (!fgets(buf, sizeof(buf)-1, attr_file)) {
    fprintf(stderr, "ERROR: could not read header line in attriubutes file %s (%s)\n",
            attr_filename, strerror(errno));
    return -1;
  }
  token = strtok_r(buf, delims, &saveptr);
  while(token) {
    attr_names = (char **)safe_realloc(attr_names, 
                                       (num_attributes + 1) * sizeof(char *));
    attr_names[num_attributes++] = safe_strdup(token);
    token = strtok_r(NULL, delims, &saveptr);
  }
  saveptr = NULL; /* reset strtok() for next line */

  /* Now that we know how many attributes there are, allocate space for values */
  attr_values = (int **)safe_malloc(num_attributes * sizeof(int **));
  for (i = 0; i < num_attributes; i++)
    attr_values[i] = (int *)safe_malloc(num_nodes * sizeof(int));

  if (!fgets(buf, sizeof(buf)-1, attr_file)) {
    fprintf(stderr, "ERROR: could not read first values line in attributes file %s (%s)\n",
            attr_filename, strerror(errno));
    return -1;
  }
  while (!feof(attr_file)) {
    thisline_values = 0;
    token = strtok_r(buf, delims, &saveptr);
    while(token) {
      if (strcasecmp(token, NA_STRING) == 0) {
        val = isBinary ? BIN_NA : CAT_NA;
      } else {
        if (sscanf(token, "%u", &val) != 1) {
          fprintf(stderr, "ERROR: bad value '%s' for node %u\n", token, nodenum);
          return -1;
        }
        if (isBinary && (val != 0 && val != 1)) {
          fprintf(stderr, "ERROR: bad value %d for binary attribute %s on node %u in attributes file %s\n",
                  val,  thisline_values < num_attributes? attr_names[thisline_values] : "UNKNOWN", 
                  nodenum, attr_filename);
          return -1;
        } else if (!isBinary && val < 0) {
          fprintf(stderr, "ERROR: bad value %d for categorical attribute %s on node %u in attributes file %s\n",
                  val,  thisline_values < num_attributes? attr_names[thisline_values] : "UNKNOWN", 
                  nodenum, attr_filename);
          return -1;
        }
      }
      if (thisline_values < num_attributes && nodenum < num_nodes) 
        attr_values[thisline_values][nodenum] = val;
      thisline_values++;
      token = strtok_r(NULL, delims, &saveptr);
    }
    if (thisline_values != num_attributes) {
      fprintf(stderr, "ERROR: %u values for node %u but expected %u in file %s\n",
              thisline_values, nodenum, num_attributes, attr_filename);
      return -1;
    }
    if (!fgets(buf, sizeof(buf)-1, attr_file)) {
      if (!feof(attr_file)) {
        fprintf(stderr, "ERROR: attempting to read attributes in file %s (%s)\n",
                attr_filename, strerror(errno));
        return -1;
      }
    }
    nodenum++;
    saveptr = NULL; /* reset strtok() for next line */
  }
  if (nodenum != num_nodes) {
    fprintf(stderr, "ERROR: %u rows after header but expected %u in file %s\n",
            nodenum, num_nodes, attr_filename);
    return -1;
  }
  fclose(attr_file);
  *out_attr_names = attr_names;
  *out_attr_values = attr_values;
  return num_attributes;
}

/*
 * Load floating point (double) attributes from file.
 * The format of the file is a header line with whitespace
 * delimited attribute names, and each subsequent line
 * the attribute values for each attribute.
 * The first line (after the header) has the values for
 * node 0, then the next line node 1, and so on.
 * 
 * Valid values are standard C library floating point format, or
 * NA for missing data (note internally this is converted to IEEE floating
 * point NaN value, so nan entered here will also be treated as missing data).
 *
 * Parameters:
 *   attr_filenname - filename of file to read
 *   num_nodes - number of nodes (must be this many values)
 *   out_attr_names - (Out) attribute names array
 *   out_attr_values - (Out) (*attr_values)[u][i] is value of attr u for node i
 * 
 * Return value:
 *   Number of attributes, or -1 on error.
 *
 * The attribute names and values arrays are allocated by 
 * this function.
 */
static int load_float_attributes(const char *attr_filename,
                                 uint_t num_nodes,
                                 char ***out_attr_names,
                                 double ***out_attr_values)
{
  const char *delims    = " \t\r\n"; /* strtok_r() delimiters  */
  uint_t nodenum        = 0;   /* node number values are for */
  uint_t num_attributes = 0;   /* number of different attributes */
  uint_t thisline_values= 0;   /* number values read this line */
  char  **attr_names   = NULL; /* array of attribute names */
  double **attr_values = NULL; /* attr_values[u][i] is value of attr u for node i */
  char *saveptr        = NULL; /* for strtok_r() */
  char *token          = NULL; /* from strtok_r() */
  char  *endptr;               /* for strtod() */
  FILE *attr_file;
  char buf[BUFSIZE];
  uint_t  i;
  double  val;

  if (!(attr_file = fopen(attr_filename, "r"))) {
    fprintf(stderr, "ERROR: could not open continuous attribute file %s (%s)\n",
            attr_filename, strerror(errno));
    return -1;
  }
  if (!fgets(buf, sizeof(buf)-1, attr_file)) {
    fprintf(stderr, "ERROR: could not read header line in continuous attriubutes file %s (%s)\n",
            attr_filename, strerror(errno));
    return -1;
  }
  token = strtok_r(buf, delims, &saveptr);
  while(token) {
    attr_names = (char **)safe_realloc(attr_names, 
                                       (num_attributes + 1) * sizeof(char *));
    attr_names[num_attributes++] = safe_strdup(token);
    token = strtok_r(NULL, delims, &saveptr);
  }
  saveptr = NULL; /* reset strtok() for next line */

  /* Now that we know how many attributes there are, allocate space for values */
  attr_values = (double **)safe_malloc(num_attributes * sizeof(double **));
  for (i = 0; i < num_attributes; i++)
    attr_values[i] = (double *)safe_malloc(num_nodes * sizeof(double));

  if (!fgets(buf, sizeof(buf)-1, attr_file)) {
    fprintf(stderr, "ERROR: could not read first values line in continuous attributes file %s (%s)\n",
            attr_filename, strerror(errno));
    return -1;
  }
  while (!feof(attr_file)) {
    thisline_values = 0;
    token = strtok_r(buf, delims, &saveptr);
    while(token) {
      if (strcasecmp(token, NA_STRING) == 0) {
        val = NAN; /* NA value for continuous is floating point NaN */
      } else {
        val = strtod(token, &endptr);
        if (*endptr != '\0') {
          fprintf(stderr, "ERROR: bad floating point value '%s' for node %u\n", token, nodenum);
          return -1;        
        }
      }
      if (thisline_values < num_attributes && nodenum < num_nodes) 
        attr_values[thisline_values][nodenum] = val;
      thisline_values++;
      token = strtok_r(NULL, delims, &saveptr);
    }
    if (thisline_values != num_attributes) {
      fprintf(stderr, "ERROR: %u values for node %u but expected %u in file %s\n",
              thisline_values, nodenum, num_attributes, attr_filename);
      return -1;
    }
    if (!fgets(buf, sizeof(buf)-1, attr_file)) {
      if (!feof(attr_file)) {
        fprintf(stderr, "ERROR: attempting to read attributes in file %s (%s)\n",
                attr_filename, strerror(errno));
        return -1;
      }
    }
    nodenum++;
    saveptr = NULL; /* reset strtok() for next line */
  }
  if (nodenum != num_nodes) {
    fprintf(stderr, "ERROR: %u rows after header but expected %u in file %s\n",
            nodenum, num_nodes, attr_filename);
    return -1;
  }
  fclose(attr_file);
  *out_attr_names = attr_names;
  *out_attr_values = attr_values;
  return num_attributes;
}


   
/*****************************************************************************
 *
 * externally visible functions
 *
 ****************************************************************************/

/*
 * Get entry for (i, j) in hashtable.
 *
 *
 * Parameters:
 *     h - hash table
 *     i - node id of source
 *     j - node id of destination
 *
 * Return value:
 *     value for key (i, j) in hashtable or 0 if none
 */
uint_t get_twopath_entry(khash_t(m64) *h, uint_t i, uint_t j)
{
  int      is_missing;
  uint64_t kiter;
  uint64_t key = ((uint64_t)i << 32) | (j & 0xffffffff);
  kiter = kh_get(m64, h, key);
  is_missing = (kiter == kh_end(h));
  if (is_missing) {
    return 0;
  }
  else {
    return kh_value(h, kiter);
  }
}



/* 
 * Return density of graph
 * 
 * Parameters:
 *   g  - digraph 
 *
 * Return value:
 *   Density of g
 */
double density(const digraph_t *g)
{
  return (double)g->num_arcs / (double)(g->num_nodes * (g->num_nodes - 1));
}

/*
 * Test if arc i -> j exists
 *
 * Parameters:
 *   g - digraph
 *   i - node to test arc from
 *   j - node to test arc to
 *
 * Return value:
 *   TRUE iff arc i->j exists
 */
bool isArc(const digraph_t *g, uint_t i, uint_t j)
{
  uint_t k;
  assert(i < g->num_nodes);
  assert(j < g->num_nodes);
  for (k = 0; k < g->outdegree[i]; k++)  {
    if (g->arclist[i][k] == j) {
      return TRUE;
    }
  }
  return FALSE;
}

/*
 * Test if an arc in either direction between i and j (i.e. either
 *  i -> j or j -> i) exists
 *
 * Parameters:
 *   g - digraph
 *   i - node to test arc from/to
 *   j - node to test arc from/to
 *
 * Return value:
 *   TRUE iff arc i->j or j->i exists
 */
bool isArcIgnoreDirection(const digraph_t *g, uint_t i, uint_t j)
{
  return isArc(g, i, j) || isArc(g, j, i);
}

/*
 * Insert arc i -> j into digraph g, WITHOUT updating allarcs flat arc list
 *
 * Parameters:
 *   g - digraph
 *   i - node to insert arc from
 *   j - node to insert arc to
 *
 * Return value:
 *   None
 */
void insertArc(digraph_t *g, uint_t i, uint_t j)
{
  assert(i < g->num_nodes);
  assert(j < g->num_nodes);
  g->num_arcs++;
  g->arclist[i] = (uint_t *)safe_realloc(g->arclist[i],
                                         (g->outdegree[i]+1) * sizeof(uint_t));
  g->arclist[i][g->outdegree[i]++] = j;
  g->revarclist[j] = (uint_t *)safe_realloc(g->revarclist[j],
                                            (g->indegree[j]+1) * sizeof(uint_t));
  g->revarclist[j][g->indegree[j]++] = i;
  updateTwoPathsMatrices(g, i, j, TRUE);
  DIGRAPH_DEBUG_PRINT(("insertArc %u -> %u indegree(%u) = %u outdegre(%u) = %u\n", i, j, j, g->indegree[j], i, g->outdegree[i]));
  assert(isArc(g, i, j));

  /* update zone information for snowball conditional estimation */
  if (g->zone[i] > g->zone[j]) {
    assert(g->zone[i] == g->zone[j] + 1);
    g->prev_wave_degree[i]++;
  } else if (g->zone[j] > g->zone[i]) {
    assert(g->zone[j] == g->zone[i] + 1);
    g->prev_wave_degree[j]++;
  }
}

/*
 * Remove arc i -> j from digraph g, WITHOUT updating allarcs flat arc list
 *
 * Parameters:
 *   g - digraph
 *   i - node to remove arc from
 *   j - node to remove arc to
 *
 * Return value:
 *   None
 */
void removeArc(digraph_t *g, uint_t i, uint_t j)
{
  uint_t k;
  DIGRAPH_DEBUG_PRINT(("removeArc %u -> %u indegree(%u) = %u outdegre(%u) = %u\n", i, j, j, g->indegree[j], i, g->outdegree[i]));
  assert(isArc(g, i, j));
  assert(i < g->num_nodes);
  assert(j < g->num_nodes);
  assert(g->num_arcs > 0);
  assert(g->outdegree[i] > 0);
  assert(g->indegree[j] > 0);
#ifdef ORDERED_ARCLIST
  /* we have to find the entry for j in the arc list and then move 
     everything after it back one, overwriting it. Note if it is
     last entry then memmove() has size 0 but that is fine.
     Same for reverse arc list. */
  for (k = 0; k < g->outdegree[i] && g->arclist[i][k] != j; k++)
    /*nothing*/;
  assert(g->arclist[i][k] == j);
  memmove(&g->arclist[i][k], &g->arclist[i][k+1],
          sizeof(uint_t)*(g->outdegree[i]-k));
  //     This is needed to keep the arclist ordered.
  //     However  this code is only used when ORDERD_ARCLIST is
  //     defined, which it is not, as I decided we did not need it. (We do
  //     not keep the arc list in any order -- but I tested the code to
  //     do it in case we do want to do so in the future).
  for (k = 0; k < g->indegree[j] && g->revarclist[j][k] != i; k++)
    /*nothing*/;
  assert(g->revarclist[j][k] == i);
  memmove(&g->revarclist[j][k], &g->revarclist[j][k+1],
          sizeof(uint_t)*(g->indegree[j]-k));
#else
  /* arclist is not ordered, so  just replace deleted entry
     with last entry */
  for (k = 0; k < g->outdegree[i] && g->arclist[i][k] != j; k++)
    /*nothing*/;
  assert(g->arclist[i][k] == j);
  g->arclist[i][k] = g->arclist[i][g->outdegree[i]-1];
  for (k = 0; k < g->indegree[j] && g->revarclist[j][k] != i; k++)
    /*nothing*/;
  assert(g->revarclist[j][k] == i);
  g->revarclist[j][k] = g->revarclist[j][g->indegree[j]-1];
#endif

  g->num_arcs--;
  g->outdegree[i]--;
  g->indegree[j]--;
  updateTwoPathsMatrices(g, i, j, FALSE);

  /* update zone information for snowball conditional estimation */ 
  if (g->zone[i] > g->zone[j]) {
    assert(g->prev_wave_degree[i] > 1);
    g->prev_wave_degree[i]--;
  } else if (g->zone[j] > g->zone[i]) {
    assert(g->prev_wave_degree[j] > 1);
    g->prev_wave_degree[j]--;
  }
}


/*
 * Insert arc i -> j into digraph g, updating allarcs flat arc list
 *
 * Parameters:
 *   g - digraph
 *   i - node to insert arc from
 *   j - node to insert arc to
 *
 * Return value:
 *   None
 */
void insertArc_allarcs(digraph_t *g, uint_t i, uint_t j)
{
  insertArc(g, i, j);
  g->allarcs = (nodepair_t *)safe_realloc(g->allarcs,
                                          g->num_arcs * sizeof(nodepair_t));
  g->allarcs[g->num_arcs-1].i = i;
  g->allarcs[g->num_arcs-1].j = j;
}

/*
 * Remove arc i -> j from digraph g, updating allarcs flat arc list
 *
 * Parameters:
 *   g - digraph
 *   i - node to remove arc from
 *   j - node to remove arc to
 *   arcidx - index in allarcs flat arc list of the i->j entry for fast removal
 *            as this is known (arc has been selected from this list)
 *
 * Return value:
 *   None
 */
void removeArc_allarcs(digraph_t *g, uint_t i, uint_t j, uint_t arcidx)
{
  removeArc(g, i, j);
  /* remove entry from the flat all arcs list */
  assert(g->allarcs[arcidx].i == i && g->allarcs[arcidx].j == j);
  /* replace deleted entry with last entry */
  /* g->num_arcs already decremented by removeArc() */
  g->allarcs[arcidx].i = g->allarcs[g->num_arcs].i;
  g->allarcs[arcidx].j = g->allarcs[g->num_arcs].j;
}


/*
 * Allocate the  digraph structure for empty digraph with given
 * number of nodes.
 *
 * Parameters:
 *    num_vertices - number of nodes in digraph
 *
 * Return values:
 *    Allocated and initizlied to empty digraph
 */
digraph_t *allocate_digraph(uint_t num_vertices)
{
  digraph_t *g = (digraph_t *)safe_malloc(sizeof(digraph_t));  
  g->num_nodes = num_vertices;
  g->num_arcs = 0;
  g->outdegree = (uint_t *)safe_calloc(num_vertices, sizeof(uint_t));
  g->arclist = (uint_t **)safe_calloc(num_vertices, sizeof(uint_t *));
  g->indegree = (uint_t *)safe_calloc(num_vertices, sizeof(uint_t));
  g->revarclist = (uint_t **)safe_calloc(num_vertices, sizeof(uint_t *));
  g->allarcs = NULL;

  MEMUSAGE_DEBUG_PRINT(("Allocated two arrays %f MB each for adjacency lists\n",
                        ((double)sizeof(uint_t) * num_vertices) / (1024*1024)));

  /* TODO change dense matrices to sparse (hash table or CSR etc.) for scalabiity */  
  g->mixTwoPathMatrix = (uint_t *)safe_calloc(num_vertices * num_vertices,
                                              sizeof(uint_t));
  g->inTwoPathMatrix = (uint_t *)safe_calloc(num_vertices * num_vertices,
                                             sizeof(uint_t));
  g->outTwoPathMatrix = (uint_t *)safe_calloc(num_vertices * num_vertices,
                                              sizeof(uint_t));

  MEMUSAGE_DEBUG_PRINT(("Allocated three arrays of %f MB each for two-paths\n",
                        ((double)sizeof(uint_t) * num_vertices * num_vertices) /
                        (1024*1024)));

  g->mixTwoPathHashTab = kh_init(m64);
  
  g->num_binattr = 0;
  g->binattr_names = NULL;
  g->binattr = NULL;
  g->num_catattr = 0;
  g->catattr_names = NULL;
  g->catattr = NULL;
  g->num_contattr = 0;
  g->contattr_names = NULL;
  g->contattr = NULL;

  g->zone  = (uint_t *)safe_calloc(num_vertices, sizeof(uint_t));
  g->max_zone = 0;
  g->num_inner_nodes = 0;
  g->inner_nodes = NULL;
  g->prev_wave_degree  = (uint_t *)safe_calloc(num_vertices, sizeof(uint_t));
  return g;
}

/*
 * Free the digraph internal structures and digraph itself
 *
 * Parameters:
 *    g - digraph to deallocate
 * Return value:
 *    None
 * Note the pointer g itelf is freed in this function
 */
void free_digraph(digraph_t *g)
{
  uint_t i;

  for (i = 0; i < g->num_binattr; i++) {
    free(g->binattr_names[i]);
    free(g->binattr[i]);
  }
  free(g->binattr);
  free(g->binattr_names);
  for (i = 0; i < g->num_catattr; i++) {
    free(g->catattr_names[i]);
    free(g->catattr[i]);
  }
  free(g->catattr);
  free(g->catattr_names);
  for (i = 0; i < g->num_contattr; i++) {
    free(g->contattr_names[i]);
    free(g->contattr[i]);
  }
  free(g->contattr);
  free(g->contattr_names);
  for (i = 0; i < g->num_nodes; i++)  {
    free(g->arclist[i]);
    free(g->revarclist[i]);
  }
  free(g->allarcs);
  free(g->arclist);
  free(g->revarclist);
  free(g->indegree);
  free(g->outdegree);
  free(g->mixTwoPathMatrix);
  free(g->inTwoPathMatrix);
  free(g->outTwoPathMatrix);
  kh_destroy(m64, g->mixTwoPathHashTab);
  free(g->zone);
  free(g->inner_nodes);
  free(g->prev_wave_degree);
  free(g);
}



/*
 * Build digraph from Pajek format arc list file,
 * with optional binary and categorical attributes files.
 *
 * In the Pajek format *vertices at top, then followed by one line for each
 * vertex (just vertex number) then *arcs followed by arcs list one per
 * line. In this program the nodes must be numbered 1..N.
 *
 * The format of the attriubtes files is header line with whitespace-delimited
 * attribute names, followed by (whitespace delimited) attributes 
 * one line per node (corresponding to node number order).
 * If the attribute file handles are NULL then no attributes.
 *
 * Parameters:
 *    pajek_file   - Pajek format arclist file handle (open read).
 *                   Closed by this function at end.
 *    binattr_filename - binary attributes filebane or NULL
 *    catattr_filename - categorical attribute filename or NULL
 *    contattr_filename- continuous attribule filename or NULL
 *
 * Return value:
 *    digraph object built from files.
 *
 * Note this function calls exit() on error.
 */
digraph_t *load_digraph_from_arclist_file(FILE *pajek_file,
                                          const char *binattr_filename,
                                          const char *catattr_filename,
                                          const char *contattr_filename)
{
  digraph_t *g = NULL;
  int i, j;
  char *p;
  char *saveptr   = NULL; /* for strtok_r() */
  char *token     = NULL; /* from strtok_r() */
  const char *delims = " \t\r\n"; /* strtok_r() delimiters for header lines */
  int num_vertices = 0;
  int num_attr;
#ifdef DEBUG_MEMUSAGE
  uint_t k, total_degree = 0;
#endif /* DEBUG_MEMUSAGE */

  char buf[BUFSIZE];
  /* the first lines should be e.g.
   * *vertices 36
   * for Pajek format
   */
  fgets(buf, sizeof(buf)-1, pajek_file);
  for (p = buf; *p !='\0'; p++) {
    *p = tolower(*p);
  }
  if (sscanf(buf, "*vertices %d\n", &num_vertices) != 1) {
    fprintf(stderr, "ERROR: expected *vertices n line but didn't find it\n");
    exit(1);
  }
  g = allocate_digraph(num_vertices);
  
  do {
    fgets(buf, sizeof(buf)-1, pajek_file);
  } while (!feof(pajek_file) && strncasecmp(buf, "*arcs", 5) != 0);
  if (feof(pajek_file)) {
    fprintf(stderr, "did not find *arcs line\n");
    exit(1);
  }
  if (!fgets(buf, sizeof(buf)-1, pajek_file)) {
    fprintf(stderr, "ERROR: attempting to read first arc  (%s)\n", strerror(errno));
    fclose(pajek_file);
    exit(1);
  }
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
    if (!isArc(g, i, j)){
      insertArc_allarcs(g, i, j); /* also update flat arclist allarcs */
    }
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
  for (k = 0; k < g->num_nodes; k++) {
    total_degree += g->outdegree[k];
  }
  MEMUSAGE_DEBUG_PRINT(("Allocated additional %f MB (twice) for %u arcs\n",
                        ((double)sizeof(uint_t) * total_degree) / (1024*1024),
                         g->num_arcs));

    MEMUSAGE_DEBUG_PRINT(("MixTwoPath hash table has %u entries (approx. %f MB) which is %f%% nonzero in dense matrix\n",
                          kh_size(g->mixTwoPathHashTab),
                          (double)(kh_size(g->mixTwoPathHashTab)*2*
                                   sizeof(uint64_t))/(1024*1024),
                          100*(double)kh_size(g->mixTwoPathHashTab) /
                          (g->num_nodes*g->num_nodes)));
#endif /* DEBUG_MEMUSAGE */
  
  if (binattr_filename) {
    if ((num_attr = load_integer_attributes(binattr_filename, num_vertices,
                                            TRUE, &g->binattr_names,
                                            &g->binattr)) < 0){
      fprintf(stderr, "ERROR: loading binary attributes from file %s failed\n", 
              binattr_filename);
      exit(1);
    }
    g->num_binattr = (uint_t)num_attr;
  }

  if (catattr_filename) {
    if ((num_attr = load_integer_attributes(catattr_filename, num_vertices,
                                            FALSE, &g->catattr_names,
                                            &g->catattr)) < 0){
      fprintf(stderr, "ERROR: loading categorical attributes from file %s failed\n", 
              catattr_filename);
      exit(1);
    }
    g->num_catattr = (uint_t)num_attr;
  }
  if (contattr_filename) {
    if ((num_attr = load_float_attributes(contattr_filename, num_vertices,
                                          &g->contattr_names,
                                          &g->contattr)) < 0){
      fprintf(stderr, "ERROR: loading continuous attributes from file %s failed\n", 
              contattr_filename);
      exit(1);
    }
    g->num_contattr = (uint_t)num_attr;
  }  
  
  
  return(g);
}

/*
 * Write arc list to stdout
 *
 * Parameters:
 *     g - digraph to dump
 *
 * Return value:
 *    None.
 *
 */
void dump_digraph_arclist(const digraph_t *g)
{
  write_digraph_arclist_to_file(stdout, g);
}

/*
 * Write some summary statistics of graph and attribute data to stdout
 */
void print_data_summary(const digraph_t * g)
{
  uint_t i,j;
  uint_t num_na_values;
  
  printf("Digraph with %u vertics and %u arcs (density %g)\n",
         g->num_nodes, g->num_arcs, density(g));
  printf("%u binary attributes\n", g->num_binattr);
  for (i = 0; i < g->num_binattr; i++) {
    printf("  %s", g->binattr_names[i]);
    num_na_values = 0;
    for (j = 0; j < g->num_nodes; j++) {
      if (g->binattr[i][j] == BIN_NA) {
        num_na_values++;
      }
    }
    printf(" has %u NA values\n", num_na_values);
  }
  printf("%u categorical attributes\n", g->num_catattr);
  for (i = 0; i < g->num_catattr; i++) {
    printf("  %s", g->catattr_names[i]);
    num_na_values = 0;
    for (j = 0; j < g->num_nodes; j++) {
      if (g->catattr[i][j] == CAT_NA) {
        num_na_values++;
      }
    }
    printf(" has %u NA values\n", num_na_values);
  }
  printf("%u continuous attributes\n", g->num_contattr);
  for (i = 0; i < g->num_contattr; i++) {
    printf("  %s", g->contattr_names[i]);
    num_na_values = 0;
    for (j = 0; j < g->num_nodes; j++) {
      if (isnan(g->contattr[i][j])) {
        num_na_values++;
      }
    }
    printf(" has %u NA values\n", num_na_values);
  }
}

/*
 * Write some statistics about the snowball sampling zones to stdout.
 */
void print_zone_summary(const digraph_t *g)
{
  uint_t   i;
  uint_t  *zone_sizes; /* number of nodes in each zone */
  uint_t   num_zones = g->max_zone + 1;

  if (num_zones == 1) {
    printf("No zone information (all nodes in zone 0)\n");
    return;
  }
  zone_sizes = (uint_t *)safe_calloc(num_zones, sizeof(uint_t));
  for (i = 0; i < g->num_nodes; i++) {
    assert(g->zone[i] < num_zones);
    zone_sizes[g->zone[i]]++;
  }
  printf("Number of zones: %u (%u waves)\n", num_zones, num_zones-1);
  printf("Number of nodes in inner waves: %u\n", g->num_inner_nodes);  
  printf("Number of nodes in each zone:\n");
  for (i = 0; i < num_zones; i++) {
    printf(" %u: %u\n", i, zone_sizes[i]);
  }
  
  free(zone_sizes);
}

/*
 * Write arc list in Pajek format to file. The node numbers are
 * 1..n in this format (not 0..n-1).
 *
 * Parameters:
 *     fp - open (write) file pointer to write to
 *     g - digraph to dump
 *
 * Return value:
 *    None.
 *
 */
void write_digraph_arclist_to_file(FILE *fp, const digraph_t *g)
{
  uint_t i, j, count=0;

  fprintf(fp, "*vertices %u\n", g->num_nodes);
  for (i = 0; i < g->num_nodes; i++)
    fprintf(fp, "%u\n", i+1);
  fprintf(fp, "*arcs\n");
  for (i = 0; i < g->num_nodes; i++)  {
    for (j = 0; j < g->outdegree[i]; j++) {
      count++;
      fprintf(fp, "%u %u\n", i+1, g->arclist[i][j]+1); /* output is 1 based */
      assert(isArc(g, i, g->arclist[i][j]));
    }
  }
  assert(count == g->num_arcs);
}


/*
 * Read snowball sampling zone file and put zone information in digraph g
 *
 * Parameters:
 *    g             - (in/out) digraph to put zone information in
 *    zone_filename - filename of zone file to read.
 *
 * Return value:
 *    0 if OK else nonzero for error.
 * 
 * The zone, max_zone, num_inner_nodes, inner_nodes, and
 * prev_wave_degree fields of g are set here.
 *
 * The format of the file is the same as that for categorical
 * attributes (and the same function is used to parse it): a header
 * line which must have just the name "zone", and each subsequent line
 * the the snowball sampling zone for each node.  The first line (after the
 * header) has the value for node 0, then the next line node 1, and
 * so on. The zones are numbered from 0 for the seed nodes.
 * 
 * E.g.:
 *
 * zone
 * 0
 * 1
 * 1
 * 2
 */
int add_snowball_zones_to_digraph(digraph_t *g, const char *zone_filename)
{
  int      num_attr, j;
  char   **attr_names;
  int    **zones;
  uint_t   i, u, v;
  uint_t  *zone_sizes; /* number of nodes in each zone */
  uint_t   num_zones;

  
  if ((num_attr = load_integer_attributes(zone_filename, g->num_nodes,
                                          FALSE, &attr_names,
                                          &zones)) < 0){
    fprintf(stderr, "ERROR: loading zones from file %s failed\n", 
            zone_filename);
    return -1;
  }
  if (num_attr != 1) {
    fprintf(stderr, "ERROR: expecting only zone attribute in zone file %s "
            "but found %d attributes\n", zone_filename, num_attr);
    return -1;
  }
  if (strcasecmp(attr_names[0], "zone") != 0) {
    fprintf(stderr, "ERROR: expecting only zone attribute in zone file %s "
            " but found %s\n", zone_filename, attr_names[0]);
    return -1;
  }
  for (i = 0; i < g->num_nodes; i++) {
    g->zone[i] = zones[0][i];
    if (g->zone[i] > g->max_zone) {
      g->max_zone = g->zone[i];
    }
  }

  num_zones = g->max_zone + 1;

  /* check that the zones are not invalid, no skipped zones */
  zone_sizes = (uint_t *)safe_calloc(num_zones, sizeof(uint_t));
  for (i = 0; i < g->num_nodes; i++) {
    assert(g->zone[i] < num_zones);
    zone_sizes[g->zone[i]]++;
  }
  for (i = 0; i < num_zones; i++) {
    if (zone_sizes[i] == 0) {
      fprintf(stderr,
              "ERROR: Max zone is %u but there are no nodes in zone %u\n",
              g->max_zone, i);
      return -1;
    }
  }

  /*
   * For conditional estimation, the zone of each node is fixed, as
   * well as all the ties between nodes in the outermost wave (last
   * zone) and ties from nodes in the last zone to nodes in the
   * second-last zone. So in MCMC procedure we to need find nodes only
   * in the inner waves (i.e. all those apart from the outermost). So
   * to all this done to be done efficiently we build the inner_nodes
   * array which is an array of size num_inner_nodes (the number of
   * nodes in zones other than the last) of each node id in an inner
   * zone.
   */
  for (i = 0; i < g->max_zone; i++) {
    g->num_inner_nodes += zone_sizes[i];
  }
  g->inner_nodes = (uint_t *)safe_calloc(g->num_inner_nodes, sizeof(uint_t));
  for (u = 0, i = 0; u < g->num_nodes; u++) {
    if (g->zone[u] < g->max_zone) {
      assert(i < g->num_inner_nodes);
      g->inner_nodes[i++] = u;
    }
  }
  
  /*
   * build prev_wave_degree[] which for each node gives the number of
   * edges to or from (i.e. ignoring direction of arc) that node
   * to/from nodes in the immediately preceding zone. (This value will always
   * be zero for all seed nodes i.e. nodes in zone 0).
   */
  for (i = 0; i < g->num_arcs; i++) {
    u = g->allarcs[i].i;
    v = g->allarcs[i].j;
    if (g->zone[u] != g->zone[v] &&
        g->zone[u] != g->zone[v] + 1 && g->zone[v] != g->zone[u] + 1){
      fprintf(stderr, "ERROR: invalid snowball zones for adjacent nodes %u "
              "(zone %u) and %u (zone %u)\n", u, g->zone[u], v, g->zone[v]);
      return -1;
    }
    if (g->zone[u] > g->zone[v]) {
      assert(g->zone[u] == g->zone[v] + 1);
      g->prev_wave_degree[u]++;
    } else if (g->zone[v] > g->zone[u]) {
      assert(g->zone[v] == g->zone[u] + 1);
      g->prev_wave_degree[v]++;
    }
  }
  
  
  for (j = 0; j < num_attr; j++) {
    free(attr_names[j]);
    free(zones[j]);
  }
  free(attr_names);
  free(zones);
  free(zone_sizes);
  return 0;
}


/*
 * Write snowball sampling zone information (used for conditional estimation)
 * to stdout for debugging.
 *
 * Parmaters:
 *   g - digraph object to dump zone info from
 *
 * Return value:
 *  None
 */
void dump_zone_info(const digraph_t *g)
{
  uint_t   i;
  uint_t   num_zones = g->max_zone + 1;

  if (num_zones == 1) {
    printf("No zone information (all nodes in zone 0)\n");
    return;
  }
  printf("Number of zones: %u (%u waves)\n", num_zones, num_zones-1);
  printf("Number of nodes in inner waves: %u\n", g->num_inner_nodes);
  printf("Nodes in inner waves:");
  for (i = 0; i < g->num_inner_nodes; i++) {
    printf(" %u", g->inner_nodes[i]);
  }
  printf("\n");
  printf("Wave of each node:");
  for (i = 0; i < g->num_nodes; i++) {
    printf(" %u", g->zone[i]);
  }
  printf("\n");
  printf("Number of ties to/from previous wave for each node:");
  for (i = 0; i < g->num_nodes; i++) {
    printf(" %u", g->prev_wave_degree[i]);
  }
  printf("\n");
}
