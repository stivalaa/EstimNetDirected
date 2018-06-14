/*****************************************************************************
 * 
 * File:    digraph.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * Directed graph data structure. Stored as arc lists (both forward and
 * a "reversed" version, for fast iteration over both in- and out- neighbours)
 * and fast lookup matrices for two-paths.
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
void updateTwoPathsMatrices(digraph_t *g, uint_t i, uint_t j, bool isAdd)
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
  }
  for (k = 0; k < g->outdegree[j]; k++) {
    v = g->arclist[j][k];
    if (v == i || v == j)
      continue;
    assert(isArc(g,j,v));
    g->mixTwoPathMatrix[INDEX2D(i, v, g->num_nodes)] += incval;
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
 * Insert arc i -> j into digraph g 
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
}

/*
 * Remove arc i -> j from digraph g 
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
  /* TODO change dense matrices to sparse (hash table or CSR etc.) for scalabiity */  
  g->mixTwoPathMatrix = (uint_t *)safe_calloc(num_vertices * num_vertices,
                                              sizeof(uint_t));
  g->inTwoPathMatrix = (uint_t *)safe_calloc(num_vertices * num_vertices,
                                             sizeof(uint_t));
  g->outTwoPathMatrix = (uint_t *)safe_calloc(num_vertices * num_vertices,
                                              sizeof(uint_t));

  g->num_binattr = 0;
  g->binattr_names = NULL;
  g->binattr = NULL;
  g->num_catattr = 0;
  g->catattr_names = NULL;
  g->catattr = NULL;
  g->num_contattr = 0;
  g->contattr_names = NULL;
  g->contattr = NULL;

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
  free(g->arclist);
  free(g->revarclist);
  free(g->indegree);
  free(g->outdegree);
  free(g->mixTwoPathMatrix);
  free(g->inTwoPathMatrix);
  free(g->outTwoPathMatrix);
  free(g);
}



/*
 * Build digraph from Pajek format arc list file,
 * with optional binary and categorical attributes files.
 *
 * In the Pajek format *vertices at top, then followed by one line for each
 * vertex (just vertex number) then *arcs followed by arcs list one per
 * line .
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
      insertArc(g, i, j);
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
  uint_t i, j, count=0;

  printf("*vertices %u\n", g->num_nodes);
  printf("*arcs\n");
  for (i = 0; i < g->num_nodes; i++)  {
    for (j = 0; j < g->outdegree[i]; j++) {
      count++;
      printf("%u %u\n", i, g->arclist[i][j]);
      assert(isArc(g, i, g->arclist[i][j]));
    }
  }
  assert(count == g->num_arcs);
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
