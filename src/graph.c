/*****************************************************************************
 * 
 * File:    digraph.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * Directed or undirected graph data structure.  Also handles bipartite
 * (two-mode)x
 *
 * For directed, stored as arc lists (both forward and a "reversed"
 * version, for fast iteration over both in- and out- neighbours).
 * Also, fast lookup hash tables for two-paths, and flat arcs or edges
 * list for fast finding of random arc.
 *
 * Preprocessor defines used:
 *
 *    TWOPATH_LOOKUP      - use two-path lookup tables (arrays by default)
 *    TWOPATH_HASHTABLES  - use hash tables (only if TWOPATH_LOOKUP defined)
 *    ORDERED_ARCLIST     - keep arclists sorted (not completely implemented)
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
#include "graph.h"


/*****************************************************************************
 *
 * local constants
 *
 ****************************************************************************/

static const size_t BUFSIZE = 16384;  /* line buffer size for reading files */
static const char *NA_STRING = "NA"; /* string in attributes files to indicate
                                        missing data (case insensitive) */
static const char *SET_NONE_STRING = "NONE"; /* string in set attribute file
                                                to indicate no elements in
                                                set (case insensitive) */

   
/*****************************************************************************
 *
 * local functions
 *
 ****************************************************************************/




#ifdef TWOPATH_LOOKUP
#ifdef TWOPATH_HASHTABLES
/*
 * Update entry for (i, j) in hashtable.
 *
 * Used for the two-path hash tables to count two-paths between nodes
 * i and j. Much more efficient to use a hash table (or some other way
 * of storing a sparse matrix) as typically at most only on the order of 1% to
 * 10% of entries are nonzero, and for large networks (on the order
 * of a million nodes) this can be orders of magnitude smaller still e.g.
 * 0.01% to 0.1%. In such cases storing in dense matrix format
 * is entirely infeasible.
 *
 * For bipartite (two-mode) networks, i is an A node and j is a B node.
 *
 * Parameters:
 *     h - pointer to hash table (pointer itself)
 *     i - node id of source
 *     j - node id of destination
 *     incval - value to add to existing value (or insert if not exists)
 *              NB this can be negative (it is either -1 or +1)
 *
 * Return value:
 *     None.
 */
static void update_twopath_entry(twopath_record_t **h, uint_t i, uint_t j,
                                 int incval)
{
  twopath_record_t rec;
  twopath_record_t *newrec;
  twopath_record_t *p;

  /* there is no reason this function cannot have other values of incval
     but in updateTwoPathsMatrices() only +1 or -1 is used */
  assert(incval == 1 || incval == -1);
  
  memset(&rec, 0, sizeof(twopath_record_t));
  rec.key.i = i;
  rec.key.j = j;
  HASH_FIND(hh, *h, &rec.key, sizeof(nodepair_t), p);
  if (p) {
    p->value += incval;
    if (p->value == 0) {
      assert(incval < 0); /* value added must have been -ve to get to zero */
      /* delete entry with value 0 to save memory (get_twopath_entry() 
         returns 0 for value if key not in hash table, so having key in table
         with value 0 is the same as having it not present).
         If we don't do this then the hash table will just keep growing
         and could use far more memory than if we do this */
      HASH_DELETE(hh, *h, p);
      free(p);
    }
  } else {
    newrec = (twopath_record_t *)safe_malloc(sizeof(*newrec));
    newrec->key.i = i;
    newrec->key.j = j;
    newrec->value = incval;
    HASH_ADD(hh, *h, key, sizeof(nodepair_t), newrec);
  }
} 
#endif /* TWOPATH_HASHTABLES */


#ifdef TWOPATH_HASHTABLES
/*
 * Update the two-paths hash tables used for fast computation of change
 * statistics for either adding or removing arc i->j
 * The hash tables in the digraph are updated in-place
 *
 * For bipartite (two-mode) networks, i is an A node and j is a B node.
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
static void updateTwoPathsMatrices(graph_t *g, uint_t i, uint_t j, bool isAdd)
{
  uint_t v,k;
  int incval = isAdd ? 1 : -1;

  if (g->is_directed) {
    assert(!g->is_bipartite); /* directed bipartite not handled yet */
    for (k = 0; k < g->outdegree[i]; k++) {
      v = g->arclist[i][k];
      if (v == i || v == j)
        continue;
      /*removed as slows significantly: assert(isArc(g,i,v)); */
      update_twopath_entry(&g->outTwoPathHashTab, v, j, incval);
      update_twopath_entry(&g->outTwoPathHashTab, j, v, incval);
    }
    for (k = 0; k < g->indegree[j]; k++) {
      v = g->revarclist[j][k];
      if (v == i || v == j)
        continue;
      /*removed as slows significantly: assert(isArc(g,v,j)); */
      update_twopath_entry(&g->inTwoPathHashTab, v, i, incval);
      update_twopath_entry(&g->inTwoPathHashTab, i, v, incval);
    }
    for (k = 0; k < g->indegree[i]; k++)  {
      v = g->revarclist[i][k];
      if (v == i || v == j)
        continue;
      /*removed as slows significantly: assert(isArc(g,v,i));*/
      update_twopath_entry(&g->mixTwoPathHashTab, v, j, incval);
    }
    for (k = 0; k < g->outdegree[j]; k++) {
      v = g->arclist[j][k];
      if (v == i || v == j)
        continue;
      /*removed as slows significantly: assert(isArc(g,j,v));*/
      update_twopath_entry(&g->mixTwoPathHashTab, i, v, incval);
    }
  } else if (g->is_bipartite) {
    /* bipartite */
    assert(!g->is_directed); /* directed bipartite not handled yet */
    assert(bipartite_node_mode(g, i) == MODE_A &&
	   bipartite_node_mode(g, j) == MODE_B);
    for (k = 0; k < g->degree[i]; k++)  {
      v = g->edgelist[i][k];
      assert(bipartite_node_mode(g, v) == MODE_B);
      if (v == j)
        continue;
      update_twopath_entry(&g->twoPathHashTabB, j, v, incval);
      update_twopath_entry(&g->twoPathHashTabB, v, j, incval);
    }
    for (k = 0; k < g->degree[j]; k++)  {
      v = g->edgelist[j][k];
      assert(bipartite_node_mode(g, v) == MODE_A);
      if (v == i)
        continue;
      update_twopath_entry(&g->twoPathHashTabA, i, v, incval);
      update_twopath_entry(&g->twoPathHashTabA, v, i, incval);
    }
  } else {
    /* undirected (one-mode) */
    for (k = 0; k < g->degree[i]; k++)  {
      v = g->edgelist[i][k];
      if (v == i || v == j)
        continue;
      update_twopath_entry(&g->twoPathHashTab, v, j, incval);
      update_twopath_entry(&g->twoPathHashTab, j, v, incval);
    }
    for (k = 0; k < g->degree[j]; k++)  {
      v = g->edgelist[j][k];
      if (v == i || v == j)
        continue;
      update_twopath_entry(&g->twoPathHashTab, v, i, incval);
      update_twopath_entry(&g->twoPathHashTab, i, v, incval);
    }
  }
}
#else /* using arrays not hash tables for two-path lookup */
/*
 * Update the two-paths matrices used for fast computation of change
 * statistics for either adding or removing arc i->j
 * The matcies in the digraph are updated in-place
 *
 * For bipartite (two-mode) networks, i is an A node and j is a B node.
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
static void updateTwoPathsMatrices(graph_t *g, uint_t i, uint_t j, bool isAdd)
{
  uint_t v,k;
  int incval = isAdd ? 1 : -1;

  if (g->is_directed) {
    assert(!g->is_directed); /* directed bipartite not handled yet */
    for (k = 0; k < g->outdegree[i]; k++) {
      v = g->arclist[i][k];
      if (v == i || v == j)
        continue;
      /*removed as slows significantly: assert(isArc(g,i,v)); */
      g->outTwoPathMatrix[INDEX2D(v, j, g->num_nodes)] += incval;
      g->outTwoPathMatrix[INDEX2D(j, v, g->num_nodes)] += incval;
    }
    for (k = 0; k < g->indegree[j]; k++) {
      v = g->revarclist[j][k];
      if (v == i || v == j)
        continue;
      /*removed as slows significantly: assert(isArc(g,v,j)); */
      g->inTwoPathMatrix[INDEX2D(v, i, g->num_nodes)] += incval;
      g->inTwoPathMatrix[INDEX2D(i, v, g->num_nodes)] += incval;
    }
    for (k = 0; k < g->indegree[i]; k++)  {
      v = g->revarclist[i][k];
      if (v == i || v == j)
        continue;
      /*removed as slows significantly: assert(isArc(g,v,i));*/
      g->mixTwoPathMatrix[INDEX2D(v, j, g->num_nodes)]+=incval;
    }
    for (k = 0; k < g->outdegree[j]; k++) {
      v = g->arclist[j][k];
      if (v == i || v == j)
        continue;
      /*removed as slows significantly: assert(isArc(g,j,v));*/
      g->mixTwoPathMatrix[INDEX2D(i, v, g->num_nodes)] += incval;
    }
  } else if (g->is_bipartite) {
    assert(!g->is_directed); /* directed bipartite not handled yet */
    assert(bipartite_node_mode(g, i) == MODE_A &&
	   bipartite_node_mode(g, j) == MODE_B);
    for (k = 0; k < g->degree[i]; k++)  {
      v = g->edgelist[i][k];
      assert(bipartite_node_mode(g, v) == MODE_B);
      if (v == j)
        continue;
      /* Note subtracting num_A_nodes as B nodes are numbered
	 num_A_nodes .. num_nodes */
      g->twoPathMatrixB[INDEX2D(j-g->num_A_nodes, v-g->num_A_nodes, g->num_nodes)]+=incval;
      g->twoPathMatrixB[INDEX2D(v-g->num_A_nodes, j-g->num_A_nodes, g->num_nodes)]+=incval;
    }
    for (k = 0; k < g->degree[j]; k++)  {
      v = g->edgelist[j][k];
      assert(bipartite_node_mode(g, v) == MODE_A);
      if (v == i)
        continue;
      g->twoPathMatrixA[INDEX2D(i, v, g->num_nodes)]+=incval;
      g->twoPathMatrixA[INDEX2D(v, i, g->num_nodes)]+=incval;
    }
  } else {
    /* undirected (one-mode) */
    for (k = 0; k < g->degree[i]; k++)  {
      v = g->edgelist[i][k];
      if (v == i || v == j)
        continue;
      g->twoPathMatrix[INDEX2D(v, j, g->num_nodes)]+=incval;
      g->twoPathMatrix[INDEX2D(j, v, g->num_nodes)]+=incval;
    }
    for (k = 0; k < g->degree[j]; k++)  {
      v = g->edgelist[j][k];
      if (v == i || v == j)
        continue;
      g->twoPathMatrix[INDEX2D(v, i, g->num_nodes)]+=incval;
      g->twoPathMatrix[INDEX2D(i, v, g->num_nodes)]+=incval;
    }
  }
}
#endif /*TWOPATH_HASHTABLES*/
#endif /*TWOPATH_LOOKUP */

#ifdef TWOPATH_LOOKUP
#ifdef TWOPATH_HASHTABLES
/*
 * Delete all entries and entire hash table.
 *
 * Parameters:
 *   h - hash table to destroy
 *
 * Return value:
 *  None.
 */
static void deleteAllHashTable(twopath_record_t *h)
{
#ifdef DO_DELETE_HASH_ENTRIES
  twopath_record_t *curr, *tmp;
  HASH_ITER(hh, h, curr, tmp) {
    HASH_DEL(h, curr);
    free(curr);
  }
#else
  (void)h; /* suppress unused parameter warning */
  /* actually deleting all the entries can take a lot of time
     and there is really no point since this is only called on exit anyway
     so do nothing here */
#endif /*DO_DELETE_HASH_ENTRIES*/
}
#endif /*TWOPATH_HASHTABLES*/
#endif /*TWOPATH_LOOKUP*/

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
    fprintf(stderr, "ERROR: could not read header line in attribuutes file %s (%s)\n",
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
  attr_values = (int **)safe_malloc(num_attributes * sizeof(int *));
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
    fprintf(stderr, "ERROR: could not read header line in continuous attributes file %s (%s)\n",
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
  attr_values = (double **)safe_malloc(num_attributes * sizeof(double *));
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


/*
 * Load set (of categorical) attributes from file.
 * The format of the file is a header line with whitespace
 * delimited attribute names, and each subsequent line
 * the attribute values for each attribute.
 * The format of the attribute values (set of categories) is a comma
 * delimited (note must have no whitespace as different attributes
 * are delimited by whitespace) list of integers making up the set of
 * categories.
 * The first line (after the header) has the values for
 * node 0, then the next line node 1, and so on.
 * 
 * E.g.:
 *
 * type   class
 * 0,9    1
 * 3,1    2,3,4,5,6
 * 2      3,10,98
 * none   0
 *
 *
 * Valid values are integer >= 0 for each of the comma-delimited categories,
 * or a single NONE (case insensitive) for empty set, 
 * or a single NA (case insensitve) for missing data.
 *
 * The highest value of any integer for an attribute gives the size of
 * the set for that attribute. The values do not need to be contiguous,
 * and the set is stored as an array of set_elem_e for maximum flexibility
 * (rather than more efficient fixed size bit set),
 * so e.g. for 'type' in the example above the set is an array of size
 * 10 indexed 0..9 as 9 is the highest value, and for 'class' an array
 * of size 99 indexed 0..98 as 98 is the highest value.
 * 
 * Note that NONE results simply in all elements of the set being
 * absent with the normal semantics of the set, however NA results in
 * all elements of the set (array) being set to SET_ELEM_NA meaning there is
 * really a single NA for that set attribute on that node, there is
 * no individual meaning of SET_ELEM_NA at a particular index in the array.
 *
 *
 * Parameters:
 *   attr_filenname - filename of file to read
 *   num_nodes - number of nodes (must be this many values)
 *   out_attr_names - (Out) attribute names array
 *   out_attr_values - (Out) (*attr_values)[u][i] is value of attr u for node i
 *   out_set_sizes   - (Out) size of set for each attribute
 * 
 * Return value:
 *   Number of attributes, or -1 on error.
 *
 * The attribute names and values arrays are allocated by 
 * this function.
 */
static int load_set_attributes(const char   *attr_filename,
                               uint_t        num_nodes,
                               char        ***out_attr_names,
                               set_elem_e ****out_attr_values,
                               uint_t       **out_set_sizes)
{
  const char *delims    = " \t\r\n"; /* strtok_r() delimiters  */
  uint_t nodenum        = 0;   /* node number values are for */
  uint_t num_attributes = 0;   /* number of different attributes */
  uint_t thisline_values= 0;   /* number values read this line */
  char  **attr_names   = NULL; /* array of attribute names */
  set_elem_e ***attr_values = NULL; /* attr_values[u][i] is value of attr u for node i */
  char *saveptr        = NULL; /* for strtok_r() */
  char *token          = NULL; /* from strtok_r() */
  uint_t  *setsizes    = NULL; /* max integer in set for each attribute */
  FILE *attr_file;
  char buf[BUFSIZE];
  uint_t  i;
  set_elem_e   *setval = NULL;
  int     pass;
  bool    firstpass;

  if (!(attr_file = fopen(attr_filename, "r"))) {
    fprintf(stderr, "ERROR: could not open set attribute file %s (%s)\n",
            attr_filename, strerror(errno));
    return -1;
  }
  if (!fgets(buf, sizeof(buf)-1, attr_file)) {
    fprintf(stderr, "ERROR: could not read header line in set attributes file %s (%s)\n",
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
  setsizes = (uint_t *)safe_calloc((size_t)num_attributes, sizeof(uint_t));
  attr_values = (set_elem_e ***)safe_malloc(num_attributes * sizeof(set_elem_e **));
  for (i = 0; i < num_attributes; i++)
    attr_values[i] = (set_elem_e **)safe_malloc(num_nodes * sizeof(set_elem_e *));

  for (pass = 0; pass < 2; pass++) {
    /* on first pass, get max int in set for each attribute so can allocate
       arrays, actually build them on second pass */
    firstpass = (pass == 0);
    if (!firstpass) {
      /* on second pass have to reopen file and skip over header line */
      if (!(attr_file = fopen(attr_filename, "r"))) {
        fprintf(stderr, "ERROR: could not open set attribute file %s (%s)\n",
            attr_filename, strerror(errno));
        return -1;
      }
      if (!fgets(buf, sizeof(buf)-1, attr_file)) {
        fprintf(stderr, "ERROR: could not read header line in set attributes file %s (%s)\n",
                attr_filename, strerror(errno));
        return -1;
      }
      GRAPH_DEBUG_PRINT(("load_set_attributes pass %d reopen %s at '%s'\n",
                           pass, attr_filename, buf));
    }
    if (!fgets(buf, sizeof(buf)-1, attr_file)) {
      fprintf(stderr, "ERROR: could not read first values line in set attributes file %s (%s)\n",
              attr_filename, strerror(errno));
      return -1;
    }
    nodenum = 0;
    saveptr = NULL; /* reset strtok() for next line */
    while (!feof(attr_file)) {
      thisline_values = 0;
      token = strtok_r(buf, delims, &saveptr);
      while(token) {
        GRAPH_DEBUG_PRINT(("load_set_attributes pass %u token '%s'\n",
                             pass, token));
        if (!firstpass) {
          setval = (set_elem_e *)safe_malloc(setsizes[thisline_values] *
                                       sizeof(set_elem_e));
        }
        if (parse_category_set(token, firstpass,
                               &setsizes[thisline_values], setval) < 0) {
          fprintf(stderr, "ERROR: bad set value '%s' for node %u\n", token,
                  nodenum);
          return -1;
        }
        if (thisline_values < num_attributes && nodenum < num_nodes) {
          if (!firstpass) {
            attr_values[thisline_values][nodenum] = setval;
          }
        }
        thisline_values++;
        token = strtok_r(NULL, delims, &saveptr);
      }
      if (thisline_values != num_attributes) {
        fprintf(stderr, "ERROR: %u set values for node %u but expected %u in file %s\n",
                thisline_values, nodenum, num_attributes, attr_filename);
        return -1;
      }
      if (!fgets(buf, sizeof(buf)-1, attr_file)) {
        if (!feof(attr_file)) {
          fprintf(stderr, "ERROR: attempting to read set attributes in file %s (%s)\n",
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
  }
  *out_attr_names = attr_names;
  *out_attr_values = attr_values;
  *out_set_sizes = setsizes;
  return num_attributes;
}

   
/*****************************************************************************
 *
 * externally visible functions
 *
 ****************************************************************************/

#ifdef TWOPATH_LOOKUP
#ifdef TWOPATH_HASHTABLES
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
uint_t get_twopath_entry(twopath_record_t *h, uint_t i, uint_t j)
{
  twopath_record_t rec, *p = NULL;
  rec.key.i = i;
  rec.key.j = j;
  HASH_FIND(hh, h, &rec.key, sizeof(nodepair_t), p);
  return (p ? p->value : 0);
}
#endif /*TWOPATH_HASHTABLES*/
#else /* not using two-path lookup tables (either arrays or hashtables) */

/* In these functions we have to count paths i -- v -- j for different
   directions (i.e -- can be <- or ->) depending on the function. So we
   can iterate over neighbours of i or of j in the outermost loop. The 
   result is the same but it can be faster if the outermost loop is over
   the node with smallest degree. Note that in the inner loop we could also
   iterate neighbours of v or of i or j, but in practice I found it far
   faster to always just use j or i rather than v for large network with very
   high maximum degree and very skewed degree distribution (physician referral
   network). On most networks (smaller, less skewed, lower max degree) it
   makes no real difference */

/* 
 * Count two-paths for (i, j): paths  i -> v -> j for some v
 */
uint_t mixTwoPaths(const graph_t *g, uint_t i, uint_t j)
{
  uint_t v,k,l;
  uint_t count = 0;

  assert(g->is_directed);
  
  if (g->outdegree[i] < g->indegree[j]) {
    for (k = 0; k < g->outdegree[i]; k++)  {
      v = g->arclist[i][k];   /* i -> v */
      if (v == i || v == j)
        continue;
      for (l = 0; l < g->indegree[j]; l++) {
        if (g->revarclist[j][l] == v) {   /* v -> j */
          count++;
        }
      }
    }
  } else {
    for (k = 0; k < g->indegree[j]; k++) {
      v = g->revarclist[j][k];  /* v -> j */
      if (v == i || v == j)
        continue;
      for (l = 0; l < g->outdegree[i]; l++) {
        if (g->arclist[i][l] == v) { /* i -> v */
          count++;
        }
      }
    }
  }
  return count;
}

/* 
 * Count out-two-paths for (i, j): paths  i <- v -> j for some v
 */
uint_t outTwoPaths(const graph_t *g, uint_t i, uint_t j)
{
  uint_t v,k,l;
  uint_t count = 0;

  assert(g->is_directed);
  
  if (g->indegree[i] < g->indegree[j]) {
    for (k = 0; k < g->indegree[i]; k++)  {
      v = g->revarclist[i][k];   /* i <- v */
      if (v == i || v == j)
        continue;
      for (l = 0; l < g->indegree[j]; l++) {
        if (g->revarclist[j][l] == v) {   /* v -> j */
          count++;
        }
      }
    }
  } else {
    for (k = 0; k < g->indegree[j]; k++) {
      v = g->revarclist[j][k]; /* v -> j */
      if (v == i || v == j) 
        continue;
      for (l = 0; l < g->indegree[i]; l++) {
        if (g->revarclist[i][l] == v) { /* i <- v */
          count++;
        }
      }
    }
  }
  return count;
}

/* 
 * Count in-two-paths for (i, j): paths  i -> v <- j for some v
 */
uint_t inTwoPaths(const graph_t *g, uint_t i, uint_t j)
{
  uint_t v,k,l;
  uint_t count = 0;

  assert(g->is_directed);
  
  if (g->outdegree[i] < g->outdegree[j]) {
    for (k = 0; k < g->outdegree[i]; k++)  {
      v = g->arclist[i][k];   /* i -> v */
      if (v == i || v == j)
        continue;
      for (l = 0; l < g->outdegree[j]; l++) {
        if (g->arclist[j][l] == v) {   /* v <- j */
          count++;
        }
      }
    }
  } else {
    for (k = 0; k < g->outdegree[j]; k++) {
      v = g->arclist[j][k]; /* v <- j */
      if (v == i || v == j)
        continue;
      for (l = 0; l < g->outdegree[i]; l++) {
        if (g->arclist[i][l] == v) { /* i -> v */
          count++;
        }
      }
    }
  }
  return count;
}


/* 
 * Count undirected two-paths for (i, j): paths  i -- v -- j for some v
 */
uint_t twoPaths(const graph_t *g, uint_t i, uint_t j)
{
  uint_t v,k,l;
  uint_t count = 0;

  assert(!g->is_directed);
  if (g->degree[i] < g->degree[j]) {
    for (k = 0; k < g->degree[i]; k++)  {
      v = g->edgelist[i][k];   /* i -- v */
      if (v == i || v == j)
        continue;
      for (l = 0; l < g->degree[j]; l++) {
        if (g->edgelist[j][l] == v) {   /* v -- j */
          count++;
        }
      }
    }
  } else {
    for (k = 0; k < g->degree[j]; k++) {
      v = g->edgelist[j][k];  /* v -- j */
      if (v == i || v == j)
        continue;
      for (l = 0; l < g->degree[i]; l++) {
        if (g->edgelist[i][l] == v) { /* i -- v */
          count++;
        }
      }
    }
  }
  return count;
}


#endif /*TWOPATH_LOOKUP*/

/*
 * Return number of potential arcs or edges in a graph
 *
 * Parameters:
 *   g          - graph or digraph
 *   allowLoops - allow self-edges (loops)
 *
 * Return value:
 *   Number of potential arcs or edges in g
 */
double num_graph_dyads(const graph_t *g, bool allowLoops)
{
  if (g->is_bipartite) {
    assert(!allowLoops);
    assert(!g->is_directed); /* directed bipartite not supported yet */
    return ((double)g->num_A_nodes * g->num_B_nodes)/2.0;
  } else {
    if (g->is_directed) {
      if (allowLoops)
	return (double)g->num_nodes * g->num_nodes;
      else
	return (double)g->num_nodes * (g->num_nodes - 1);
    } else {
      if (allowLoops)
	return ((double)g->num_nodes * g->num_nodes)/2.0;
      else
	return ((double)g->num_nodes * (g->num_nodes - 1))/2.0;
    }
  }
}

/*
 * Return number of potential arcs or edges between inner snowball nodes a graph
 * Must not allow loops
 *
 * Parameters:
 *   g          - graph or digraph
 *
 * Return value:
 *   Number of potential arcs or edges between inner nodes g
 */

double num_graph_inner_dyads(const graph_t *g)
{
  assert(!g->is_bipartite); /* snowball sampled bipartite not supported */
  return g->is_directed ?
    (double)g->num_inner_nodes*(g->num_inner_nodes-1) :
    (double)g->num_inner_nodes*(g->num_inner_nodes-1)/2.0;
}

/*
 * Return number of arcs (digraph) or edges (graph)
 *
 * Parameters:
 *   g          - graph or digraph
 *
 * Return value:
 *   Number of edges or arcs
 */
uint_t num_arcs_or_edges(const graph_t *g)
{
  if (g->is_directed)
    return g->num_arcs;
  else
    return g->num_edges;
}

/*
 * Return number of inner snowball zone arcs (digraph) or edges (graph)
 *
 * Parameters:
 *   g          - graph or digraph
 *
 * Return value:
 *   Number of inner snowball zone edges or arcs
 */
uint_t num_inner_arcs_or_edges(const graph_t *g)
{
  assert(!g->is_bipartite); /* snowball sampled bipartite not supported */
  if (g->is_directed)
    return g->num_inner_arcs;
  else
    return g->num_inner_edges;
}

/*
 * Return density of graph
 * 
 * Parameters:
 *   g          - graph or digraph
 *   allowLoops - allow self-edges (loops)
 *
 * Return value:
 *   Density of g
 */
double density(const graph_t *g, bool allowLoops)
{
  return (double)num_arcs_or_edges(g) / num_graph_dyads(g, allowLoops);
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
bool isArc(const graph_t *g, uint_t i, uint_t j)
{
  uint_t k;
  assert(g->is_directed);
  assert(i < g->num_nodes);
  assert(j < g->num_nodes);
  if (g->outdegree[i] < g->indegree[j]) {
    for (k = 0; k < g->outdegree[i]; k++)  {
      if (g->arclist[i][k] == j) {
        return TRUE;
      }
    }
  } else {
    for (k = 0; k < g->indegree[j]; k++) {
      if (g->revarclist[j][k] == i) {
        return TRUE;
      }
    }
  }
  return FALSE;
}

/*
 * Test if edge i -- j exists
 *
 * Parameters:
 *   g - graph
 *   i - node to test edge
 *   j - node to test edge
 *
 * Return value:
 *   TRUE iff edge i -- j exists
 */
bool isEdge(const graph_t *g, uint_t i, uint_t j)
{
  uint_t k;
  assert(!g->is_directed);
  assert(i < g->num_nodes);
  assert(j < g->num_nodes);
  if (g->degree[i] < g->degree[j]) {
    for (k = 0; k < g->degree[i]; k++)  {
      if (g->edgelist[i][k] == j) {
        return TRUE;
      }
    }
  } else {
    for (k = 0; k < g->degree[j]; k++) {
      if (g->edgelist[j][k] == i) {
        return TRUE;
      }
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
bool isArcIgnoreDirection(const graph_t *g, uint_t i, uint_t j)
{
  return g->is_directed ? (isArc(g, i, j) || isArc(g, j, i)) : isEdge(g, i, j);
}

/*
 * Test if an arc i -> j in digraph or edge i -- j in undirected graph exists
 *
 * Parameters:
 *   g - graph or digraph
 *   i - node to test arc from or edge
 *   j - node to test arc to or edge
 *
 * Return value:
 *   TRUE iff : for digraph arc i->j exists, 
 *              for graph edge i--j exists.
 */
bool isArcOrEdge(const graph_t *g, uint_t i, uint_t j)
{
  return g->is_directed ? isArc(g, i, j) : isEdge(g, i, j);
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
void insertArc(graph_t *g, uint_t i, uint_t j)
{
  assert(g->is_directed);
  assert(i < g->num_nodes);
  assert(j < g->num_nodes);
  g->num_arcs++;
  g->arclist[i] = (uint_t *)safe_realloc(g->arclist[i],
                                         (g->outdegree[i]+1) * sizeof(uint_t));
  g->arclist[i][g->outdegree[i]++] = j;
  g->revarclist[j] = (uint_t *)safe_realloc(g->revarclist[j],
                                            (g->indegree[j]+1) * sizeof(uint_t));
  g->revarclist[j][g->indegree[j]++] = i;
#ifdef TWOPATH_LOOKUP
  updateTwoPathsMatrices(g, i, j, TRUE);
#endif /* TWOPATH_LOOKUP */
  GRAPH_DEBUG_PRINT(("insertArc %u -> %u indegree(%u) = %u outdegre(%u) = %u\n", i, j, j, g->indegree[j], i, g->outdegree[i]));
  /*removed as slows significantly: assert(isArc(g, i, j));*/

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
 * Insert edge i -- j into graph g, WITHOUT updating alledges flat edge list
 *
 * Parameters:
 *   g - graph
 *   i - node to insert edge
 *   j - node to insert edge
 *
 * Return value:
 *   None
 */
void insertEdge(graph_t *g, uint_t i, uint_t j)
{
  assert(!g->is_directed);
  assert(i < g->num_nodes);
  assert(j < g->num_nodes);
  if (g->is_bipartite) {
    assert(bipartite_node_mode(g, i) != bipartite_node_mode(g, j));
  }
  g->num_edges++;
  g->edgelist[i] = (uint_t *)safe_realloc(g->edgelist[i],
                                         (g->degree[i]+1) * sizeof(uint_t));
  g->edgelist[i][g->degree[i]++] = j;
  g->edgelist[j] = (uint_t *)safe_realloc(g->edgelist[j],
                                            (g->degree[j]+1) * sizeof(uint_t));
  g->edgelist[j][g->degree[j]++] = i;
#ifdef TWOPATH_LOOKUP
  updateTwoPathsMatrices(g, i, j, TRUE);
#endif /* TWOPATH_LOOKUP */
  GRAPH_DEBUG_PRINT(("insertEdge %u -- %u degree(%u) = %u degree(%u) = %u\n", i, j, j, g->degree[j], i, g->degree[i]));

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
void removeArc(graph_t *g, uint_t i, uint_t j)
{
  uint_t k;
  assert(g->is_directed);
  GRAPH_DEBUG_PRINT(("removeArc %u -> %u indegree(%u) = %u outdegre(%u) = %u\n", i, j, j, g->indegree[j], i, g->outdegree[i]));
  /*removed as slows significantly: assert(isArc(g, i, j));*/
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
#ifdef TWOPATH_LOOKUP
  updateTwoPathsMatrices(g, i, j, FALSE);
#endif /* TWOPATH_LOOKUP */

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
 * Remove edge i -- j from graph g, WITHOUT updating alledges flat edge list
 *
 * Parameters:
 *   g - graph
 *   i - node to remove arc from
 *   j - node to remove arc to
 *
 * Return value:
 *   None
 */
void removeEdge(graph_t *g, uint_t i, uint_t j)
{
  uint_t k;
  assert(!g->is_directed);
  GRAPH_DEBUG_PRINT(("removeEdge %u -> %u degree(%u) = %u degree(%u) = %u\n", i, j, j, g->degree[j], i, g->degree[i]));
  assert(i < g->num_nodes);
  assert(j < g->num_nodes);
  assert(g->num_edges > 0);
  assert(g->degree[i] > 0);
  assert(g->degree[j] > 0);

  /* edgelist is not ordered, so  just replace deleted entry
     with last entry */
  for (k = 0; k < g->degree[i] && g->edgelist[i][k] != j; k++)
    /*nothing*/;
  assert(g->edgelist[i][k] == j);
  g->edgelist[i][k] = g->edgelist[i][g->degree[i]-1];
  for (k = 0; k < g->degree[j] && g->edgelist[j][k] != i; k++)
    /*nothing*/;
  assert(g->edgelist[j][k] == i);
  g->edgelist[j][k] = g->edgelist[j][g->degree[j]-1];

  g->num_edges--;
  g->degree[i]--;
  g->degree[j]--;
#ifdef TWOPATH_LOOKUP
  updateTwoPathsMatrices(g, i, j, FALSE);
#endif /* TWOPATH_LOOKUP */

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
void insertArc_allarcs(graph_t *g, uint_t i, uint_t j)
{
  assert(g->is_directed);
  insertArc(g, i, j);
  g->allarcs = (nodepair_t *)safe_realloc(g->allarcs,
                                          g->num_arcs * sizeof(nodepair_t));
  g->allarcs[g->num_arcs-1].i = i;
  g->allarcs[g->num_arcs-1].j = j;
}

/*
 * Insert edge i -- j into graph g, updating alledeges flat edge list
 *
 * Parameters:
 *   g - graph
 *   i - node to insert edge
 *   j - node to insert edge
 *
 * Return value:
 *   None
 */
void insertEdge_alledges(graph_t *g, uint_t i, uint_t j)
{
  assert(!g->is_directed);
  insertEdge(g, i, j);
  g->alledges = (nodepair_t *)safe_realloc(g->alledges,
                                          g->num_edges * sizeof(nodepair_t));
  g->alledges[g->num_edges-1].i = i;
  g->alledges[g->num_edges-1].j = j;
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
void removeArc_allarcs(graph_t *g, uint_t i, uint_t j, uint_t arcidx)
{
  assert(g->is_directed);
  removeArc(g, i, j);
  /* remove entry from the flat all arcs list */
  assert(g->allarcs[arcidx].i == i && g->allarcs[arcidx].j == j);
  /* replace deleted entry with last entry */
  /* g->num_arcs already decremented by removeArc() */
  g->allarcs[arcidx].i = g->allarcs[g->num_arcs].i;
  g->allarcs[arcidx].j = g->allarcs[g->num_arcs].j;
}

/*
 * Remove edge i -- j from graph g, updating alledges flat edge list
 *
 * Parameters:
 *   g - graph
 *   i - node to remove edge
 *   j - node to remove edge
 *   edgeidx - index in aledges flat edge list of the i--j entry for
 *            fast removal
 *            as this is known (edge has been selected from this list)
 *
 * Return value:
 *   None
 */
void removeEdge_alledges(graph_t *g, uint_t i, uint_t j, uint_t edgeidx)
{
  assert(!g->is_directed);
  removeEdge(g, i, j);
  /* remove entry from the flat all edges list */
  assert(g->alledges[edgeidx].i == i && g->alledges[edgeidx].j == j);
  /* replace deleted entry with last entry */
  /* g->num_edges already decremented by removeEdge() */
  g->alledges[edgeidx].i = g->alledges[g->num_edges].i;
  g->alledges[edgeidx].j = g->alledges[g->num_edges].j;
}


/*
 * Insert arc i -> j into digraph g, updating allinnerarcs flat arc list
 *
 * Used for conditional estimation when we must add an arc that is
 * between nodes in inner zones and in same zone or adjacent zones only.
 *
 * Parameters:
 *   g - digraph
 *   i - node to insert arc from
 *   j - node to insert arc to
 *
 * Return value:
 *   None
 */
void insertArc_allinnerarcs(graph_t *g, uint_t i, uint_t j)
{
  assert(g->is_directed);
  assert(g->zone[i] < g->max_zone && g->zone[j] < g->max_zone);
  assert(labs((long)g->zone[i] - (long)g->zone[j]) <= 1);
  insertArc(g, i, j);
  g->num_inner_arcs++;
  g->allinnerarcs = (nodepair_t *)safe_realloc(g->allinnerarcs,
                                               g->num_inner_arcs *
                                               sizeof(nodepair_t));
  g->allinnerarcs[g->num_inner_arcs-1].i = i;
  g->allinnerarcs[g->num_inner_arcs-1].j = j;
}

/*
 * Insert edge i -- j into graph g, updating allinneredges flat edge list
 *
 * Used for conditional estimation when we must add an edge that is
 * between nodes in inner zones and in same zone or adjacent zones only.
 *
 * Parameters:
 *   g - graph
 *   i - node to insert edge
 *   j - node to insert edge
 *
 * Return value:
 *   None
 */
void insertEdge_allinneredges(graph_t *g, uint_t i, uint_t j)
{
  assert(!g->is_directed);
  assert(g->zone[i] < g->max_zone && g->zone[j] < g->max_zone);
  assert(labs((long)g->zone[i] - (long)g->zone[j]) <= 1);
  insertEdge(g, i, j);
  g->num_inner_edges++;
  g->allinneredges = (nodepair_t *)safe_realloc(g->allinneredges,
                                               g->num_inner_edges *
                                               sizeof(nodepair_t));
  g->allinneredges[g->num_inner_edges-1].i = i;
  g->allinneredges[g->num_inner_edges-1].j = j;
}

/*
 * Remove arc i -> j from digraph g, updating allinnerarcs flat arc list
 *
 * Used for conditional estimation when we must delete an arc that is
 * between nodes in inner zones and in same zone or adjacent zones only.
 *
 * Parameters:
 *   g - digraph
 *   i - node to remove arc from
 *   j - node to remove arc to
 *   arcidx - index in allinnerarcs flat arc list of the i->j entry for fast 
 *            removal as this is known (arc has been selected from this list)
 *
 * Return value:
 *   None
 */
void removeArc_allinnerarcs(graph_t *g, uint_t i, uint_t j, uint_t arcidx)
{
  assert(g->is_directed);
  assert(g->zone[i] < g->max_zone && g->zone[j] < g->max_zone);
  assert(labs((long)g->zone[i] - (long)g->zone[j]) <= 1);
  removeArc(g, i, j);
  /* remove entry from the flat all arcs list */
  assert(g->allinnerarcs[arcidx].i == i && g->allinnerarcs[arcidx].j == j);
  /* replace deleted entry with last entry */
  g->num_inner_arcs--;
  g->allinnerarcs[arcidx].i = g->allinnerarcs[g->num_inner_arcs].i;
  g->allinnerarcs[arcidx].j = g->allinnerarcs[g->num_inner_arcs].j;
}

/*
 * Remove edge i -- j from graph g, updating allinneredges flat edge list
 *
 * Used for conditional estimation when we must delete an edge that is
 * between nodes in inner zones and in same zone or adjacent zones only.
 *
 * Parameters:
 *   g - graph
 *   i - node to remove edge
 *   j - node to remove edge
 *   edgeidx - index in allinneredges flat edge list of the i->j entry for fast 
 *            removal as this is known (edge has been selected from this list)
 *
 * Return value:
 *   None
 */
void removeEdge_allinneredges(graph_t *g, uint_t i, uint_t j, uint_t edgeidx)
{
  assert(!g->is_directed);
  assert(g->zone[i] < g->max_zone && g->zone[j] < g->max_zone);
  assert(labs((long)g->zone[i] - (long)g->zone[j]) <= 1);
  removeEdge(g, i, j);
  /* remove entry from the flat all edges list */
  assert(g->allinneredges[edgeidx].i == i && g->allinneredges[edgeidx].j == j);
  /* replace deleted entry with last entry */
  g->num_inner_edges--;
  g->allinneredges[edgeidx].i = g->allinneredges[g->num_inner_edges].i;
  g->allinneredges[edgeidx].j = g->allinneredges[g->num_inner_edges].j;
}


/*
 * Insert arc i -> j into digraph g, updating all_maxtermsender_arcs flat arc list
 *
 * Used for citation ERGM estimation when we must add an arc that is
 * sent from a node with highest term value
 *
 * Parameters:
 *   g - digraph
 *   i - node to insert arc from
 *   j - node to insert arc to
 *
 * Return value:
 *   None
 */
void insertArc_all_maxtermsender_arcs(graph_t *g, uint_t i, uint_t j)
{
  assert(g->is_directed);
  assert(g->term[i] == g->max_term && g->term[j] <= g->max_term);
  insertArc(g, i, j);
  g->num_maxtermsender_arcs++;
  g->all_maxtermsender_arcs = (nodepair_t *)safe_realloc(g->all_maxtermsender_arcs,
                                               g->num_maxtermsender_arcs *
                                               sizeof(nodepair_t));
  g->all_maxtermsender_arcs[g->num_maxtermsender_arcs-1].i = i;
  g->all_maxtermsender_arcs[g->num_maxtermsender_arcs-1].j = j;
}

/*
 * Remove arc i -> j from digraph g, updating all_maxtermsender_arcs flat arc list
 *
 * Used for citation ERGM estimation when we must delete an arc that is
 * sent from a node with highest term value
 *
 * Parameters:
 *   g - digraph
 *   i - node to remove arc from
 *   j - node to remove arc to
 *   arcidx - index in all_maxtermsender_arcs flat arc list of the i->j entry for fast 
 *            removal as this is known (arc has been selected from this list)
 *
 * Return value:
 *   None
 */
void removeArc_all_maxtermsender_arcs(graph_t *g, uint_t i, uint_t j, uint_t arcidx)
{
  assert(g->is_directed);
  assert(g->term[i] == g->max_term && g->term[j] <= g->max_term);
  removeArc(g, i, j);
  /* remove entry from the flat all arcs list */
  assert(g->all_maxtermsender_arcs[arcidx].i == i && g->all_maxtermsender_arcs[arcidx].j == j);
  /* replace deleted entry with last entry */
  g->num_maxtermsender_arcs--;
  g->all_maxtermsender_arcs[arcidx].i = g->all_maxtermsender_arcs[g->num_maxtermsender_arcs].i;
  g->all_maxtermsender_arcs[arcidx].j = g->all_maxtermsender_arcs[g->num_maxtermsender_arcs].j;
}


/*
 * Allocate the  graph/digraph structure for empty graph/digraph with given
 * number of nodes.
 *
 * Parameters:
 *    num_vertices - number of nodes in graph/digraph
 *    is_directed  - TRUE for directed graph else undirected
 *
 * Return values:
 *    Allocated and initizlied to empty graph or digraph
 */
graph_t *allocate_graph(uint_t num_vertices, bool is_directed)
{
  graph_t *g = (graph_t *)safe_calloc(1, sizeof(graph_t));
  g->is_directed = is_directed;
  g->num_nodes = num_vertices;
  if (is_directed) {
    g->num_arcs = 0;
    g->outdegree = (uint_t *)safe_calloc((size_t)num_vertices, sizeof(uint_t));
    g->arclist = (uint_t **)safe_calloc((size_t)num_vertices, sizeof(uint_t *));
    g->indegree = (uint_t *)safe_calloc((size_t)num_vertices, sizeof(uint_t));
    g->revarclist = (uint_t **)safe_calloc((size_t)num_vertices, sizeof(uint_t *));
    g->allarcs = NULL;
    
#ifdef TWOPATH_LOOKUP
#ifdef TWOPATH_HASHTABLES
    g->mixTwoPathHashTab = NULL;
    g->inTwoPathHashTab = NULL;
    g->outTwoPathHashTab = NULL;
    
    assert(sizeof(void *) == 8); /* require 64 bit addressing for large uthash */
#ifdef DEBUG_MEMUSAGE
#ifdef HASH_BLOOM
    /* https://troydhanson.github.io/uthash/userguide.html#_bloom_filter_faster_misses */
    MEMUSAGE_DEBUG_PRINT(("Bloom filter n = %u overhead %f MB (three times)\n",
                          HASH_BLOOM, (pow(2, HASH_BLOOM)/8192)/(1024)));
#endif /* HASH_BLOOM */
#endif /* DEBUG_MEMUSAGE */
#else
    g->mixTwoPathMatrix = (uint_t *)safe_calloc((size_t)num_vertices * num_vertices,
                                                sizeof(uint_t));
    g->inTwoPathMatrix = (uint_t *)safe_calloc((size_t)num_vertices * num_vertices,
                                               sizeof(uint_t));
    g->outTwoPathMatrix = (uint_t *)safe_calloc((size_t)num_vertices * num_vertices,
                                                sizeof(uint_t));
#ifdef DEBUG_MEMUSAGE
    MEMUSAGE_DEBUG_PRINT(("mixTwoPathMatrix size %f MB\n", 
                          (double)num_vertices*num_vertices*sizeof(uint_t)/
                          (1024*1024)));
    MEMUSAGE_DEBUG_PRINT(("inTwoPathMatrix size %f MB\n", 
                          (double)num_vertices*num_vertices*sizeof(uint_t)/
                          (1024*1024)));
    MEMUSAGE_DEBUG_PRINT(("outTwoPathMatrix size %f MB\n", 
                          (double)num_vertices*num_vertices*sizeof(uint_t)/
                          (1024*1024)));
#endif /*DEBUG_MEMUSAGE*/
#endif /* TWOPATH_HASHTABLES */
#endif /* TWOPATH_LOOKUP */
  } else {
    /* undirected */
    g->num_edges = 0;
    g->degree = (uint_t *)safe_calloc((size_t)num_vertices, sizeof(uint_t));
    g->edgelist = (uint_t **)safe_calloc((size_t)num_vertices, sizeof(uint_t *));
    g->alledges = NULL;
    
#ifdef TWOPATH_LOOKUP
#ifdef TWOPATH_HASHTABLES
    g->twoPathHashTab = NULL;
    
    assert(sizeof(void *) == 8); /* require 64 bit addressing for large uthash */
#ifdef DEBUG_MEMUSAGE
#ifdef HASH_BLOOM
    /* https://troydhanson.github.io/uthash/userguide.html#_bloom_filter_faster_misses */
    MEMUSAGE_DEBUG_PRINT(("Bloom filter n = %u overhead %f MB (three times)\n",
                          HASH_BLOOM, (pow(2, HASH_BLOOM)/8192)/(1024)));
#endif /* HASH_BLOOM */
#endif /* DEBUG_MEMUSAGE */
#else
    g->twoPathMatrix = (uint_t *)safe_calloc((size_t)num_vertices * num_vertices,
                                                sizeof(uint_t));
#ifdef DEBUG_MEMUSAGE
    MEMUSAGE_DEBUG_PRINT(("twoPathMatrix size %f MB\n", 
                          (double)num_vertices*num_vertices*sizeof(uint_t)/
                          (1024*1024)));
#endif /*DEBUG_MEMUSAGE*/
#endif /* TWOPATH_HASHTABLES */
#endif /* TWOPATH_LOOKUP */
  }
  
  g->num_binattr = 0;
  g->binattr_names = NULL;
  g->binattr = NULL;
  g->num_catattr = 0;
  g->catattr_names = NULL;
  g->catattr = NULL;
  g->num_contattr = 0;
  g->contattr_names = NULL;
  g->contattr = NULL;
  g->num_setattr = 0;
  g->setattr_names = NULL;
  g->setattr_lengths = NULL;
  g->setattr = NULL;

  g->zone  = (uint_t *)safe_calloc((size_t)num_vertices, sizeof(uint_t));
  g->max_zone = 0;
  g->num_inner_nodes = 0;
  g->inner_nodes = NULL;
  g->prev_wave_degree  = (uint_t *)safe_calloc((size_t)num_vertices, sizeof(uint_t));
  g->num_inner_arcs = 0;
  g->allinnerarcs = NULL;
  g->term  = (uint_t *)safe_calloc((size_t)num_vertices, sizeof(uint_t));  
  g->max_term = 0;
  g->num_maxterm_nodes = 0;
  g->maxterm_nodes = NULL;
  g->num_maxtermsender_arcs = 0;
  g->all_maxtermsender_arcs = NULL;
  return g;
}

/*
 * Free the graph internal structures and graph itself
 *
 * Parameters:
 *    g - graph to deallocate
 * Return value:
 *    None
 * Note the pointer g itelf is freed in this function
 */
void free_graph(graph_t *g)
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
  for (i = 0; i < g->num_setattr; i++) {
    free(g->setattr_names[i]);
    free(g->setattr[i]);
  }
  free(g->setattr);
  free(g->setattr_names);
  if (g->is_directed) {
    for (i = 0; i < g->num_nodes; i++)  {
      free(g->arclist[i]);
      free(g->revarclist[i]);
    }
  } else {
    for (i = 0; i < g->num_nodes; i++)  {
      free(g->edgelist[i]);
    }
  }
  free(g->allarcs);
  free(g->arclist);
  free(g->revarclist);
  free(g->indegree);
  free(g->outdegree);
  free(g->alledges);
  free(g->edgelist);
  free(g->degree);
#ifdef TWOPATH_LOOKUP
#ifdef TWOPATH_HASHTABLES
  deleteAllHashTable(g->twoPathHashTab);
#else /* using arrays not hash tables for two-path lookup */
  free(g->twoPathMatrix);
#endif /* TWOPATH_HASHTABLES */
#endif /* TWOPATH_LOOKUP */
  free(g->zone);
  free(g->inner_nodes);
  free(g->prev_wave_degree);
  free(g->allinnerarcs);
  free(g->all_maxtermsender_arcs);
  free(g);
}


/*
 * Get number of nodes from Pajek network file.
 *
 * In the Pajek format *vertices at top, then followed by one line for each
 * vertex (just vertex number) then *arcs followed by arcs list one per
 * line. In this program the nodes must be numbered 1..N.
 *
 * Parameters:
 *    pajek_file   - Pajek format arclist file handle (open read).
 *                   Closed by this function at end.
 *
 * Return value:
 *    number of vertices as read from Pajek file.
 *
 * Note this function calls exit() on error.
 */
uint_t get_num_vertices_from_arclist_file(FILE *pajek_file)
{
  char *p;
  int num_vertices = 0;

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
  if (num_vertices < 1) {
    fprintf(stderr, "ERROR: number of vertices is %d\n", num_vertices);
    exit(1);
  }
  fclose(pajek_file);
  return (uint_t)num_vertices;
}



/*
 * Write arc list or edge list to stdout
 *
 * Parameters:
 *     g - digraph or graph to dump
 *
 * Return value:
 *    None.
 *
 */
void dump_graph_arclist(const graph_t *g)
{
  write_graph_arclist_to_file(stdout, g);
}

/*
 * Write some summary statistics of graph and attribute data to stdout
 */
void print_data_summary(const graph_t * g, bool allowLoops)
{
  uint_t i,j;
  uint_t num_na_values;
  
  printf("%s with %u vertices and %u %s (density %g) [%s]\n",
         g->is_directed ? "Digraph" : "Graph",
         g->num_nodes,
         g->is_directed ? g->num_arcs : g->num_edges,
         g->is_directed ? "arcs" : "edges",
         density(g, allowLoops),
         allowLoops ? "loops allowed" : "loops not allowed");
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
  printf("%u set attributes\n", g->num_setattr);
  for (i = 0; i < g->num_setattr; i++) {
    printf("  %s (size %u)", g->setattr_names[i], g->setattr_lengths[i]);
    num_na_values = 0;
    for (j = 0; j < g->num_nodes; j++) {
      if (g->setattr[i][j][0] == SET_ELEM_NA) {
        num_na_values++;
      }
    }
    printf(" has %u NA values\n", num_na_values);
  }
}

/*
 * Write some statistics about the snowball sampling zones to stdout.
 */
void print_zone_summary(const graph_t *g)
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
  printf("Number of arcs in inner waves: %u\n", g->num_inner_arcs);
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
void write_graph_arclist_to_file(FILE *fp, const graph_t *g)
{
  uint_t i, j, count=0;

  fprintf(fp, "*vertices %u\n", g->num_nodes);
  for (i = 0; i < g->num_nodes; i++)
    fprintf(fp, "%u\n", i+1);
  if (g->is_directed) {
    fprintf(fp, "*arcs\n");
    for (i = 0; i < g->num_nodes; i++)  {
      for (j = 0; j < g->outdegree[i]; j++) {
        count++;
        fprintf(fp, "%u %u\n", i+1, g->arclist[i][j]+1); /* output is 1 based */
        /*removed as slows significantly: assert(isArc(g, i, g->arclist[i][j]));*/
      }
    }
    assert(count == g->num_arcs);
  } else {
    /* undirected */
    fprintf(fp, "*edges\n");
    for (i = 0; i < g->num_nodes; i++)  {
      for (j = 0; j < g->degree[i]; j++) {
	if (i <= g->edgelist[i][j]) {
	  count++;
	  fprintf(fp, "%u %u\n", i+1, g->edgelist[i][j]+1); /* output is 1 based */
	}
      }
    }
    assert(count == g->num_edges);
  }
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
 * The zone, max_zone, num_inner_nodes, inner_nodes, prev_wave_degree,
 * num_inner_arcs and allinnerarcs fields of g are set here.
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
int add_snowball_zones_to_graph(graph_t *g, const char *zone_filename)
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
   * Also build allinnerarcs flat arcs list of arcs between nodes in inner waves
   * used for conditional estimation fast lookup of such an arc to delete.
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
    if (g->zone[u] < g->max_zone && g->zone[v] < g->max_zone) {
      g->num_inner_arcs++;
      GRAPH_DEBUG_PRINT(("inner arc %u: %u -> %u (zones %u %u)\n", g->num_inner_arcs-1, u, v, g->zone[i], g->zone[v]));
      g->allinnerarcs = (nodepair_t *)safe_realloc(g->allinnerarcs,
                                                   g->num_inner_arcs *
                                                   sizeof(nodepair_t));
      g->allinnerarcs[g->num_inner_arcs-1].i = u;
      g->allinnerarcs[g->num_inner_arcs-1].j = v;
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
void dump_zone_info(const graph_t *g)
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


/*
 * Parse comma-delimited list of int into set.
 *
 * The format of the input isa a comma delimited (note must have no
 * whitespace as different attributes are delimited by whitespace)
 * list of integers making up the set of categories.  Valid values are
 * integer >= 0 for each of the comma-delimited categories, or a
 * single NONE (case insensitive) for empty set, or a single NA (case
 * insensitve) for missing data.
 *
 * The highest value of any integer for an attribute gives the size of
 * the set for that attribute. The values do not need to be contiguous,
 * and the set is stored as an array of set_elem_e for maximum flexibility
 * (rather than more efficient fixed size bit set),
 *
 * Note that NONE results simply in all elements of the set being
 * absent with the normal semantics of the set, however NA results in
 * all elements of the set (array) being set to SET_ELEM_NA meaning there is
 * really a single NA for that set attribute on that node, there is
 * no individual meaning of SET_ELEM_NA at a particular index in the array.
 *
 *
 * Parameters:
 *    str       - input string comma-delimited list of nonnegative integers
 *    firstpass - if True, do not do any allocation or build output value,
 *                just set output parameter size to max int value in list
 *    size      - (in/Out) max int value parsed in set, if firstpass
 *                         else if firstpass not True then (in) size from
 *                         first pass
 *    setval    - (Out) set value parsed, if firstpass False
 *                Must be allocated by caller to large enough (from max
 *                size ever set here on firstpass). Not use if firstpass True.
 *
 * Return value:
 *    0 if OK, -1 on error.
 *
 *
 * This function has external linkage only so it can be used in unit tests.
 */
int parse_category_set(char *str, bool firstpass, uint_t *size,
                              set_elem_e *setval)
{
  const char *delims   = ","; /* strtok_r() delimiters  */
  char *saveptr        = NULL; /* for strtok_r() */
  char *token          = NULL; /* from strtok_r() */
  uint_t val            = 0;
  uint_t i;

  if (strcasecmp(str, NA_STRING) == 0) {
    if (!firstpass) {
      for (i = 0; i < *size; i++) {
        setval[i] = SET_ELEM_NA;
      }
    }
  } else if (strcasecmp(str, SET_NONE_STRING) == 0) {
    if (!firstpass) {
      for (i = 0; i < *size; i++) {
        setval[i] = SET_ELEM_ABSENT;
      }
    }
  } else {
    if (!firstpass) {
      /* on second pass set have size so all the elements to absent to start */
      for (i = 0; i < *size; i++) {
        setval[i] = SET_ELEM_ABSENT;
      }
    }
    token = strtok_r(str, delims, &saveptr);
    while(token) {
      if (sscanf(token, "%u", &val) != 1) {
        fprintf(stderr, "ERROR: bad value '%s' in set\n", token);
        return -1;
      }
      if (firstpass) {
        /* first pass, just get largest int for size of set */
        if (val + 1 > *size) {
          *size = val + 1;
          GRAPH_DEBUG_PRINT(("parse_category_set token '%s' size now %u\n",
                               token, *size));
        }
      } else {
        /* second pass, set the flags in the set for each int present in list */
        assert(val < *size);
        setval[val] = SET_ELEM_PRESENT;
      }
      token = strtok_r(NULL, delims, &saveptr);
    }
  }
  return 0;
}


/*
 * Load the nodal attributes from files.
 *
 * The format of the attributes files is header line with whitespace-delimited
 * attribute names, followed by (whitespace delimited) attributes 
 * one line per node (corresponding to node number order).
 * If the attribute file handles are NULL then no attributes.
 *
 * Parameters:
 *    g                - (in/out) digraph object
 *                       updated with attributes loaded
 *                       (number of nodes must match files)
 *    binattr_filename - binary attributes filebane or NULL
 *    catattr_filename - categorical attribute filename or NULL
 *    contattr_filename- continuous attribute filename or NULL
 *    setattr_filename - set attribute filename or NULL
 *
 * Return value:
 *    nonzero on error
 *
 */
int load_attributes(graph_t *g, 
                    const char *binattr_filename,
                    const char *catattr_filename,
                    const char *contattr_filename,
                    const char *setattr_filename)
{
  int  num_attr;
  int  i;
  bool setFailed = FALSE;
    
  if (binattr_filename) {
    if ((num_attr = load_integer_attributes(binattr_filename, g->num_nodes,
                                            TRUE, &g->binattr_names,
                                            &g->binattr)) < 0){
      fprintf(stderr, "ERROR: loading binary attributes from file %s failed\n", 
              binattr_filename);
      return 1;
    }
    g->num_binattr = (uint_t)num_attr;
  }

  if (catattr_filename) {
    if ((num_attr = load_integer_attributes(catattr_filename, g->num_nodes,
                                            FALSE, &g->catattr_names,
                                            &g->catattr)) < 0){
      fprintf(stderr, "ERROR: loading categorical attributes from file %s failed\n", 
              catattr_filename);
      return 1;
    }
    g->num_catattr = (uint_t)num_attr;
  }
  if (contattr_filename) {
    if ((num_attr = load_float_attributes(contattr_filename, g->num_nodes,
                                          &g->contattr_names,
                                          &g->contattr)) < 0){
      fprintf(stderr, "ERROR: loading continuous attributes from file %s failed\n", 
              contattr_filename);
      return 1;
    }
    g->num_contattr = (uint_t)num_attr;
  }  
  if (setattr_filename) {
    if ((num_attr = load_set_attributes(setattr_filename, g->num_nodes,
                                        &g->setattr_names,
                                        &g->setattr,
                                        &g->setattr_lengths)) < 0){
      fprintf(stderr, "ERROR: loading set attributes from file %s failed\n", 
              setattr_filename);
      return 1;
    }
    g->num_setattr = (uint_t)num_attr;
    setFailed = FALSE;
    for (i = 0; i < num_attr; i++) {
      if (g->setattr_lengths[i] == 0) {
        fprintf(stderr,
                "ERROR: all values for set attribute %s are NA or NONE\n",
                g->setattr_names[i]);
        setFailed = TRUE;
      }
      if (setFailed) {
        fprintf(stderr, "ERROR: some set attribute(s) not valid\n");
        return 1;
      }
    }
  }
  return 0;
}


/*
 * Read citation ERGM (cERGM) term file and put term information in digraph g
 *
 * Parameters:
 *    g             - (in/out) digraph to put zone information in
 *    term_filename - filename of term (time period) file to read.
 *
 * Return value:
 *    0 if OK else nonzero for error.
 *
 * The term, max_term, num_maxterm_nodes, maxterm_nodes
 * fields of g are set here.
 *
 * The format of the file is the same as that for categorical
 * attributes (and the same function is used to parse it), however
 * the integers here are interpreted as ordinal, rather than categorical.
 *
 * A header line which must have just the name "term", and each
 * subsequent line the term (time period) for each node.  The first
 * line (after the header) has the term for node 0, then the next line
 * node 1, and so on. The terms are numbered from 0 for the first term
 * (time period):
 *
 * E.g.:
 *
 * term
 * 0
 * 1
 * 1
 * 2
 */
int add_cergm_terms_to_digraph(graph_t *g, const char *term_filename)
{
  int      num_attr, k;
  char   **attr_names;
  int    **terms;
  uint_t   i,u,v;
  uint_t  *term_sizes; /* number of nodes in each term */
  uint_t   num_terms;


  if ((num_attr = load_integer_attributes(term_filename, g->num_nodes,
                                          FALSE, &attr_names,
                                          &terms)) < 0){
    fprintf(stderr, "ERROR: loading cERGM terms from file %s failed\n",
            term_filename);
    return -1;
  }
  if (num_attr != 1) {
    fprintf(stderr, "ERROR: expecting only term attribute in term file %s "
            "but found %d attributes\n", term_filename, num_attr);
    return -1;
  }
  if (strcasecmp(attr_names[0], "term") != 0) {
    fprintf(stderr, "ERROR: expecting only term attribute in term file %s "
            " but found %s\n", term_filename, attr_names[0]);
    return -1;
  }
  for (i = 0; i < g->num_nodes; i++) {
    g->term[i] = terms[0][i];
    if (g->term[i] > g->max_term) {
      g->max_term = g->term[i];
    }
  }

  num_terms = g->max_term + 1;

  /* check that the terms are not invalid, no skipped terms */
  term_sizes = (uint_t *)safe_calloc(num_terms, sizeof(uint_t));
  for (i = 0; i < g->num_nodes; i++) {
    assert(g->term[i] < num_terms);
    term_sizes[g->term[i]]++;
  }
  for (i = 0; i < num_terms; i++) {
    if (term_sizes[i] == 0) {
      fprintf(stderr,
              "ERROR: Max term is %u but there are no nodes in term %u\n",
              g->max_term, i);
      return -1;
    }
  }

  /*
   * For conditional estimation of the citation ERGM (cERGM) model,
   * the term of each node is fixed, as well as all the ties sent from
   * nodes in terms (time periods) earlier than the last (latest, most
   * recent). Only ties sent from nodes in the last (term == max_term)
   * time period can be added or deleted.
   *
   * So to do all this efficiently we build the maxterm_nodes array
   * which is an array of size num_maxterm_nodes (the number of nodes
   * in the last term) of each node id in an the last term. Then in
   * the MCMC procedure, ties can only be added or deleted if they are
   * directed from a node in this list (to any other node, in the list
   * or not).
   */
  g->num_maxterm_nodes = term_sizes[g->max_term];
  g->maxterm_nodes = (uint_t *)safe_calloc(g->num_maxterm_nodes, sizeof(uint_t));
  for (u = 0, i = 0; u < g->num_nodes; u++) {
    assert(g->term[u] <= g->max_term);
    if (g->term[u] == g->max_term) {
      assert(i < g->num_maxterm_nodes);
      g->maxterm_nodes[i++] = u;
    }
  }

  /* Similarly to snowball conditional estimation, for TNT and IFD samplers
   * build all_maxtermsender_arcs flat list of arcs sent from
   * node with maximum term value here for fast lookup of such an arc to 
   * delete.
   */
  for (i = 0; i < g->num_arcs; i++) {
    u = g->allarcs[i].i;
    v = g->allarcs[i].j;
    if (g->term[u] == g->max_term) {
      g->num_maxtermsender_arcs++;
      GRAPH_DEBUG_PRINT(("maxterm arc %u: %u -> %u (terms %u %u)\n",
                           g->num_inner_arcs-1, u, v, g->term[i], g->term[v]));
      g->all_maxtermsender_arcs =
        (nodepair_t *)safe_realloc(g->all_maxtermsender_arcs,
                                   g->num_maxtermsender_arcs *
                                   sizeof(nodepair_t));
      g->all_maxtermsender_arcs[g->num_maxtermsender_arcs-1].i = u;
      g->all_maxtermsender_arcs[g->num_maxtermsender_arcs-1].j = v;
    }
  }
  
  for (k = 0; k < num_attr; k++) {
    free(attr_names[k]);
    free(terms[k]);
  }
  free(attr_names);
  free(terms);
  free(term_sizes);
  return 0;
}


/*
 * Write cERGM term (time period)information (used for cERGM
 * conditional estimation) to stdout for debugging.
 *
 * Parmaters:
 *   g - digraph object to dump term info from
 *
 * Return value:
 *  None
 */
void dump_term_info(const graph_t *g)
{
  uint_t   i;
  uint_t   num_terms = g->max_term + 1;

  if (num_terms == 1) {
    printf("No cERGM term information (all nodes in term 0)\n");
    return;
  }
  printf("Number of cERGM terms: %u (max term %u) \n", num_terms, g->max_term);
  printf("Number of nodes in last term: %u\n", g->num_maxterm_nodes);
  printf("Nodes in last term:");
  for (i = 0; i < g->num_maxterm_nodes; i++) {
    printf(" %u", g->maxterm_nodes[i]);
  }
  printf("\n");
  printf("Term of each node:");
  for (i = 0; i < g->num_nodes; i++) {
    printf(" %u", g->term[i]);
  }
  printf("\n");
}

/*
 * Write some statistics about the cERGM time periods (terms) to stdout.
 */
void print_term_summary(const graph_t *g)
{
  uint_t   i;
  uint_t  *term_sizes; /* number of nodes in each term */
  uint_t   num_terms = g->max_term + 1;

  if (num_terms == 1) {
    printf("No cERGM term information (all nodes in term 0)\n");
    return;
  }
  term_sizes = (uint_t *)safe_calloc(num_terms, sizeof(uint_t));
  for (i = 0; i < g->num_nodes; i++) {
    assert(g->term[i] < num_terms);
    term_sizes[g->term[i]]++;
  }
  printf("Number of cERGM terms: %u (max term %u)\n", num_terms, num_terms-1);
  printf("Number of nodes in last term: %u\n", g->num_maxterm_nodes);
  printf("Number of nodes in each term:\n");
  for (i = 0; i < num_terms; i++) {
    printf(" %u: %u\n", i, term_sizes[i]);
  }
  free(term_sizes);
}

/*
 * count the number of self-edges (loops) in the graph.
 * Note: only used for warning message at start, this is not efficient,
 * especially if the change statistics for loop is being computed anyway.
 *
 * Parameters:
 *     g - digraph
 *
 * Return value:
 *     number of self-edges in g
 */
uint_t num_loops(const graph_t *g)
{
  uint_t k;
  uint_t count = 0;

  if (g->is_directed) {
    for (k = 0; k < g->num_arcs; k++) {
      if (g->allarcs[k].i == g->allarcs[k].j) {
        ++count;
      }
    }
  } else {
    for (k = 0; k < g->num_edges; k++) {
      if (g->alledges[k].i == g->alledges[k].j) {
        ++count;
      }
    }
  }
  return count;
}

/*
 * Return true if a node has a self-edge (loop)
 *
 * Parameters:
 *    g - digraph
 *    u - node in g
 *
 * Return value:
 *    True if u has a self-edge else False
 */
bool has_loop(const graph_t *g, uint_t u)
{
  uint_t k;

  if (g->is_directed) {
    for (k = 0; k < g->outdegree[u]; k++) {
      if (g->arclist[u][k] == u) {
        return TRUE;
      }
    }
  } else {
    for (k = 0; k < g->degree[u]; k++) {
      if (g->edgelist[u][k] == u) {
        return TRUE;
      }
    }
  }
  return FALSE;
}


/*
 * Insert arc i -> j into g if directed, or edge i -- j if g is undirected
 *  WITHOUT updating allarcs flat arc list
 *
 * Parameters:
 *   g - graph
 *   i - node to insert arc from
 *   j - node to insert arc to
 *
 * Return value:
 *   None
 */
void insertArcOrEdge(graph_t *g, uint_t i, uint_t j)
{
  if (g->is_directed)
    insertArc(g, i, j);
  else
    insertEdge(g, i, j);
}

/*
 * Remove arc i -> j from g if directed, or edge i -- j if g is undirected
 *  WITHOUT updating allarcs flat arc list
 *
 * Parameters:
 *   g - graph
 *   i - node to remove arc from
 *   j - node to remove arc to
 *
 * Return value:
 *   None
 */
void removeArcOrEdge(graph_t *g, uint_t i, uint_t j)
{
  if (g->is_directed)
    removeArc(g, i, j);
  else
    removeEdge(g, i, j);
}


/*
 * Insert arc i -> j into g if directed, or edge i -- j if g is undirected
 *  updating flat arc/edge list
 *
 * Parameters:
 *   g - graph
 *   i - node to insert arc from
 *   j - node to insert arc to
 *
 * Return value:
 *   None
 */
void insertArcOrEdge_updatelist(graph_t *g, uint_t i, uint_t j)
{
  if (g->is_directed)
    insertArc_allarcs(g, i, j);
  else
    insertEdge_alledges(g, i, j);
}

/*
 * Remove arc i -> j from g if directed, or edge i -- j if g is undirected
 *   updating allarcs flat arc list
 *
 * Parameters:
 *   g - graph
 *   i - node to remove arc from
 *   j - node to remove arc t
 *   idx - index in flat arc or edge list of the i->j entry for fast removal
 *            as this is known (arc/edge has been selected from this list)

 *
 * Return value:
 *   None
 */
void removeArcOrEdge_updatelist(graph_t *g, uint_t i, uint_t j, uint_t idx)
{
  if (g->is_directed)
    removeArc_allarcs(g, i, j, idx);
  else
    removeEdge_alledges(g, i, j, idx);
}



/*
 * Insert arc i -> j into g if directed, or edge i -- j if g is undirected
 *  updating inner  arc/edge list
 *
 * Used for conditional estimation when we must add an arc that is
 * between nodes in inner zones and in same zone or adjacent zones only.
 *
 * Parameters:
 *   g - graph
 *   i - node to insert arc from
 *   j - node to insert arc to
 *
 * Return value:
 *   None
 */
void insertArcOrEdge_updateinnerlist(graph_t *g, uint_t i, uint_t j)
{
  if (g->is_directed)
    insertArc_allinnerarcs(g, i, j);
  else
    insertEdge_allinneredges(g, i, j);
}

/*
 * Remove arc i -> j from g if directed, or edge i -- j if g is undirected
 *   updating allarcs flat arc list
 *
 * Used for conditional estimation when we must add an arc that is
 * between nodes in inner zones and in same zone or adjacent zones only.
 *
 * Parameters:
 *   g - graph
 *   i - node to remove arc from
 *   j - node to remove arc t
 *   idx - index in inner arc or edge list of the i->j entry for fast removal
 *            as this is known (arc/edge has been selected from this list)
 *
 * Return value:
 *   None
 */
void removeArcOrEdge_updateinnerlist(graph_t *g, uint_t i, uint_t j, uint_t idx)
{
  if (g->is_directed)
    removeArc_allinnerarcs(g, i, j, idx);
  else
    removeEdge_allinneredges(g, i, j, idx);
}


/*
 * For a bipartite (two-mode) graph, return node type (mode).
 *
 * Parameters:
 *   g - graph
 *   i - node to test
 *
 * Return value:
 *   MODE_A for node i in mode A, NODE_B for mode B
 *
 */
bipartite_node_mode_e bipartite_node_mode(const graph_t *g, uint_t i)
{
  assert(g->is_bipartite);
  assert(i < g->num_nodes);
  return (i < g->num_A_nodes ? MODE_A : MODE_B);
}
