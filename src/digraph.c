/*****************************************************************************
 * 
 * File:    digraph.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * Directed graph data structure. Stored as arc lists (both forward
 * and a "reversed" version, for fast iteration over both in- and out-
 * neighbours) and fast lookup hash tables for two-paths, and also
 * flat arcs list for fast selection of an arc uniformly at random.
 *
 *
 ****************************************************************************/

/*
 * TODO: may be better to use e.g. CSR sparse matrix format for
 * two-path fast lookup instead of hash table as we often iterate over
 * neighbours of a node i.e. along a row in computing change
 * statistics. (The same applies to the actual graph itself, currently
 * stored in adjacency lists). However then it is
 * difficult/inefficient to update when adding/deleting an arc which
 * is at the core of the MCMC algorithm so no good for that. (Note
 * also that if we stored the actual graph adjacency matrix in sparse
 * matrix format (CSR for example) then the numbers of two-paths can
 * be very efficiently computed my multiplying the matrix by itself
 * with sparse BLAS routine DCSRMM.)
 */

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
static const char *SET_NONE_STRING = "NONE"; /* string in set attribute file
                                                to indicate no elements in
                                                set (case insensitive) */

   
/*****************************************************************************
 *
 * local functions
 *
 ****************************************************************************/

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
}
#else
/*
 * Update the two-paths matrices used for fast computation of change
 * statistics for either adding or removing arc i->j
 * The matcies in the digraph are updated in-place
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
}
#endif /*TWOPATH_HASHTABLES*/

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
      DIGRAPH_DEBUG_PRINT(("load_set_attributes pass %d reopen %s at '%s'\n",
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
        DIGRAPH_DEBUG_PRINT(("load_set_attributes pass %u token '%s'\n",
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
void insertArc_allinnerarcs(digraph_t *g, uint_t i, uint_t j)
{
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
void removeArc_allinnerarcs(digraph_t *g, uint_t i, uint_t j, uint_t arcidx)
{
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
  g->outdegree = (uint_t *)safe_calloc((size_t)num_vertices, sizeof(uint_t));
  g->arclist = (uint_t **)safe_calloc((size_t)num_vertices, sizeof(uint_t *));
  g->indegree = (uint_t *)safe_calloc((size_t)num_vertices, sizeof(uint_t));
  g->revarclist = (uint_t **)safe_calloc((size_t)num_vertices, sizeof(uint_t *));
  g->allarcs = NULL;

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
  for (i = 0; i < g->num_setattr; i++) {
    free(g->setattr_names[i]);
    free(g->setattr[i]);
  }
  free(g->setattr);
  free(g->setattr_names);
  for (i = 0; i < g->num_nodes; i++)  {
    free(g->arclist[i]);
    free(g->revarclist[i]);
  }
  free(g->allarcs);
  free(g->arclist);
  free(g->revarclist);
  free(g->indegree);
  free(g->outdegree);
#ifdef TWOPATH_HASHTABLES
  deleteAllHashTable(g->mixTwoPathHashTab);
  deleteAllHashTable(g->inTwoPathHashTab);
  deleteAllHashTable(g->outTwoPathHashTab);
#else
  free(g->mixTwoPathMatrix);
  free(g->inTwoPathMatrix);
  free(g->outTwoPathMatrix);
#endif /* TWOPATH_HASHTABLES */
  free(g->zone);
  free(g->inner_nodes);
  free(g->prev_wave_degree);
  free(g->allinnerarcs);
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
 * Build digraph from Pajek format arc list file.
 * Note this builds the digraph structure only, not the node attributes.
 * These are read and built separately by calling load_attributes().
 * The graph g must be already allocated by allocate_digraph() with the
 * correct number of nodes, which can be found by
 * get_num_nodes(pajek_file).
 *
 * In the Pajek format *vertices at top, then followed by one line for each
 * vertex (just vertex number) then *arcs followed by arcs list one per
 * line. In this program the nodes must be numbered 1..N.
 *
 * Parameters:
 *    pajek_file   - Pajek format arclist file handle (open read).
 *                   Closed by this function at end.
 *   g             - (in/out) digraph object already allocated as above.
 *
 * Return value:
 *    digraph object built from files (same as parameter g)
 *
 * Note this function calls exit() on error.
 */
digraph_t *load_digraph_from_arclist_file(FILE *pajek_file, digraph_t *g)
{
  int i, j;
  char *p;
  char *saveptr   = NULL; /* for strtok_r() */
  char *token     = NULL; /* from strtok_r() */
  const char *delims = " \t\r\n"; /* strtok_r() delimiters for header lines */
  int num_vertices = 0;
#ifdef DEBUG_MEMUSAGE
  uint_t k, total_degree = 0;
#ifdef TWOPATH_HASHTABLES
  struct timeval start_timeval, end_timeval, elapsed_timeval;
  int            etime;
#endif /*TWOPATH_HASHABLES */
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

  if ((uint_t)num_vertices != g->num_nodes) {
    fprintf(stderr, "ERROR: expected %u vertices but found %d\n",
            g->num_nodes, num_vertices);
    exit(1);
  }
  
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
    /* FIXME calling isArc all the time is inefficient: since we are using
       hash table anyway would be better to check these in temporary hash
       table structure (or build graph as hash table then convert to the
       adjacency list structures). */
    if (!isArc(g, i, j)){
      insertArc_allarcs(g, i, j); /* also update flat arclist allarcs */
    }

#ifdef TWOPATH_HASHTABLES
#ifdef DEBUG_MEMUSAGE
    if (g->num_arcs % 1000 == 0){
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
#endif /* DEBUG_MEMUSAGE */

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
  
  printf("Digraph with %u vertices and %u arcs (density %g)\n",
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
      /*removed as slows significantly: assert(isArc(g, i, g->arclist[i][j]));*/
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
   * Also build allinnerarcs flat arcs list of arcs between nodes in inner waves
   * used for conditional estimation fast lookup of such an arc to delete
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
      DIGRAPH_DEBUG_PRINT(("inner arc %u: %u -> %u (zones %u %u)\n", g->num_inner_arcs-1, u, v, g->zone[i], g->zone[v]));
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
          DIGRAPH_DEBUG_PRINT(("parse_category_set token '%s' size now %u\n",
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
int load_attributes(digraph_t *g, 
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


