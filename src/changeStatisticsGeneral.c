/*****************************************************************************
 * 
 * File:    changeStatisticsGeneral.c
 * Author:  Alex Stivala
 * Created: January 2022
 *
 * Functions to compute graph change statistics, that apply to both
 * directed and undirected graphs. Each
 * function takes a pointer to a digraph struct, and two node numbers
 * i and j and returns the value of the change statistic for adding
 * the edge i -- j or arc i -> j (which must not already exist in the graph).
 *
 * Also takes lambda (decay) parameter which is only used for
 * some statistics ("alternating" statistics).
 *
 * For change statistics dependent on a nodal attribute, there is
 * an additional parameter a which is the index of the attribute
 * to use.
 *
 * On some functions there is additionally a parameter indicating when
 * the change statistic is being computed as part of a delete (rather
 * than add) move, which can be used for some implementations that can
 * be more easily implemented with this information. However in
 * general it is more elegant and simpler to compute the statistic for
 * adding the arc (for delete moves the value returned is just
 * negated, and the change statistic function does not depend on or
 * need to use this information at all).
 *
 * Some of these functions are adapted from the original PNet code by Peng Wang:
 *
 *   Wang P, Robins G, Pattison P. PNet: A program for the simulation and
 *   estimation of exponential random graph models. University of
 *   Melbourne. 2006.
 *
 * And for the definitions of the change statistics:
 * 
 *   Robins, G., Pattison, P., & Wang, P. (2009). Closure, connectivity and
 *   degree distributions: Exponential random graph (p*) models for
 *   directed social networks. Social Networks, 31(2), 105-117.
 * 
 *   Snijders, T. A., Pattison, P. E., Robins, G. L., & Handcock,
 *   M. S. (2006). New specifications for exponential random graph
 *   models. Sociological methodology, 36(1), 99-153.
 * 
 * And also generally:
 * 
 *   Lusher, D., Koskinen, J., & Robins, G. (Eds.). (2013). Exponential
 *   random graph models for social networks: Theory, methods, and
 *   applications. New York, NY: Cambridge University Press.
 * 
 * especially Ch. 6:
 *
 *   Koskinen, J., & Daraganova, G. (2013). Exponential random graph model
 *   fundamentals. In "Exponential random graph models for social networks:
 *   Theory, methods, and applications." (pp. 49-76). New York, NY:
 *   Cambridge University Press.
 *
 * As well as the statnet ergm terms, and references for specific
 * change statistics where indicated.
 *
 *
 * Do NOT compile with -ffast-math on gcc as we depend on IEEE handling of NaN
 *
 ****************************************************************************/

#include <math.h>
#include <assert.h>
#include "changeStatisticsGeneral.h"
#include "changeStatisticsBipartiteUndirected.h"


/*****************************************************************************
 *
 * local functions
 *
 ****************************************************************************/


/*
 * Size of the intersection of two sets.  Each set is represented by
 * array of set_elem_e indicating NA, ABSENT or PRESENT for each
 * element, size of intersection is number of array indices where both
 * are PRESENT. Note we ignore NA here.
 *
 * Parameters:
 *      a - set
 *      b - set
 *      n - size of arrays a and b
 *
 * Return value:
 *      number of elements that are in both sets a and b.
 */
static uint_t set_intersection_size(set_elem_e a[], set_elem_e b[], uint_t n)
{
  uint_t i;
  uint_t count = 0;

  for (i = 0; i < n; i++) {
    if (a[i] == SET_ELEM_PRESENT && b[i] == SET_ELEM_PRESENT) {
      count++;
    }
  }
  return count;
}

/*
 * Size of the union of two sets.  Each set is represented by
 * array of set_elem_e indicating NA, ABSENT or PRESENT for each
 * element, size of union is number of array indices where either (or both)
 * are PRESENT. Note we ignore NA here.
 *
 * Parameters:
 *      a - set
 *      b - set
 *      n - size of arrays a and b
 *
 * Return value:
 *      number of elements that are in either (or both) of the sets a and b.
 */
static uint_t set_union_size(set_elem_e a[], set_elem_e b[], uint_t n)
{
  uint_t i;
  uint_t count = 0;

  for (i = 0; i < n; i++) {
    if (a[i] == SET_ELEM_PRESENT || b[i] == SET_ELEM_PRESENT) {
      count++;
    }
  }
  return count;
}

/*
 * Jaccard index (similarity) for two sets. The Jaccard index is the
 * size of the intersection over the size of the union. (If a and b 
 * are both empty it is defined as 1).
 *
 * Parameters:
 *      a - set
 *      b - set
 *      n - size of arrays a and b
 *
 * Return value:
 *      Jaccard coefficient (similarity) of the two sets a and b
 *                                              
 * Has external linkage for use in unit tests.
 */
double jaccard_index(set_elem_e a[], set_elem_e b[], uint_t n)
{
  uint_t intersection_size, union_size;

  intersection_size = set_intersection_size(a, b, n);
  union_size = set_union_size(a, b, n);
  return (union_size == 0) ? 1 : (double)intersection_size / (double)union_size;
}


/*****************************************************************************
 *
 * change statistics functions
 *
 ****************************************************************************/


/************************* Structural ****************************************/

/*
 * Change statistic for Isolates
 */
double changeIsolates(graph_t *g, uint_t i, uint_t j, double lambda)
{
  double delta = 0;
  (void)lambda; /* unused parameter */
  slow_assert(!isArcOrEdge(g, i, j));
  if (g->is_directed) {
    if (g->indegree[i] == 0 && g->outdegree[i] == 0) {
      delta--;
    }
    if (i != j && g->indegree[j] == 0 && g->outdegree[j] == 0) {
      delta--;
    }
  } else {
    /* undirected */
    if (g->degree[i] == 0) {
      delta--;
    }
    if (i != j && g->degree[j] == 0) {
      delta--;
    }
  }
  return delta;
}

/*
 * Change statistic for two-path (triad census 021C; but note that since
 * these statistics, unlike motifs, are not induced subgraphs, this also
 * counts, in some cases multiple times, 111D, 111U, 030T, 030C, 201,
 * 120D, 120U, 120C, 210, 300)
 * also known as TwoMixStar (or m2star) for directed
 * and 2-star for undirected
 */
double changeTwoPath(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)lambda; /* unused parameter */
  slow_assert(!isArcOrEdge(g, i, j));

  if (g->is_directed) {
    if (i == j) {
      return 0;
    } else {
      return g->indegree[i] + g->outdegree[j] - (isArc(g, j, i) ? 2 : 0);
    }
  } else {
    /* undirected */
    if (i == j) {
      return 0;
    } else {
      return g->degree[i] + g->degree[j];
    }
  }
}


/*
 * Change statistic for loop (self-edge)
 * Note allowLoops = True must be set for this to work
 */
double changeLoop(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)g;      /* unused parameter*/
  (void)lambda; /* unused parameter */
  slow_assert(!isArcOrEdge(g, i, j));

  return i == j;
}

/************************* Actor attribute (binary) **************************/


/*
 * Change statistic for Interaction
 */
double changeInteraction(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  (void)isDelete; /* unused parameter */
  (void)exponent; /* unused parameter */
  return g->binattr[a][i] != BIN_NA && g->binattr[a][j] != BIN_NA &&
         g->binattr[a][i] && g->binattr[a][j];
}

/********************* Actor attribute (categorical) *************************/


/*
 * Change statistic for categorical matching
 */
double changeMatching(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  (void)isDelete; /* unused parameter */
  (void)exponent; /* unused parameter */
  slow_assert(!isArcOrEdge(g, i, j));
  return g->catattr[a][i] != CAT_NA && g->catattr[a][j] != CAT_NA &&
         g->catattr[a][i] == g->catattr[a][j];
}

/*
 * Change statistic for categorical mismatching
 */
double changeMismatching(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  (void)isDelete; /* unused parameter */
  (void)exponent; /* unused parameter */
  slow_assert(!isArcOrEdge(g, i, j));
  return g->catattr[a][i] != CAT_NA && g->catattr[a][j] != CAT_NA &&
         g->catattr[a][i] != g->catattr[a][j];
}


/********************* Actor attribute (continuous) *************************/

/*
 * Change statistic for continuous diff (absolute difference of attributes)
 */
double changeDiff(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  (void)isDelete; /* unused parameter */
  (void)exponent; /* unused parameter */
  slow_assert(!isArcOrEdge(g, i, j));
  if (isnan(g->contattr[a][i]) || isnan(g->contattr[a][j]))
    return 0;
  else
    return fabs(g->contattr[a][i] - g->contattr[a][j]);
}


/*
 * Change statistic for continuous sum (sum of attributes)
 */
double changeSum(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  (void)isDelete; /* unused parameter */
  (void)exponent; /* unused parameter */
  slow_assert(!isArcOrEdge(g, i, j));
  if (isnan(g->contattr[a][i]) || isnan(g->contattr[a][j]))
    return 0;
  else
    return g->contattr[a][i] + g->contattr[a][j];
}



/***************** Actor attribute (set of categorical) ********************/


/*
 * Change statistic for set Jaccard similarity
 */
double changeJaccardSimilarity(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  (void)isDelete; /* unused parameter */
  (void)exponent; /* unused parameter */
  slow_assert(!isArcOrEdge(g, i, j));
  /* For NA values all elements of set are set to NA so just check first */
  if (g->setattr[a][i][0] == SET_ELEM_NA || g->setattr[a][j][0] == SET_ELEM_NA)
    return 0;
  else
    return jaccard_index(g->setattr[a][i], g->setattr[a][j],
                         g->setattr_lengths[a]);
}


/********************* Dyadic covariate (continuous) *************************/

/* 
 * Change statistic for geographical distance between two nodes,
 * using for each node the pair of continuous attributes labelled 
 * as being latitude and longitude
 */
double changeGeoDistance(graph_t *g, uint_t i, uint_t j)
{
  double dist, lati, longi, latj, longj;
  slow_assert(!isArcOrEdge(g, i, j));

  lati  = g->contattr[g->latitude_index][i];
  longi = g->contattr[g->longitude_index][i];
  latj  = g->contattr[g->latitude_index][j];
  longj = g->contattr[g->longitude_index][j];

  if (isnan(lati) || isnan(longi) || isnan(latj) || isnan(longj)) {
    return 0;
  }
  else {
    dist = geo_distance(lati, longi, latj, longj);
    return dist;
  }
}

/* 
 * Change statistic for logarithm of geographical distance between two nodes,
 * using for each node the pair of continuous attributes labelled 
 * as being latitude and longitude
 */
double changeLogGeoDistance(graph_t *g, uint_t i, uint_t j)
{
  double dist;
  slow_assert(!isArcOrEdge(g, i, j));

  dist = changeGeoDistance(g, i, j);
  if (dist > 0) {
    return log(dist);
  } else {
    return 0;
  }
}


/* 
 * Change statistic for Euclidean distance between two nodes,
 * using for each node the triple of continuous attributes labelled 
 * as being x, y, and z coordinates.
 */
double changeEuclideanDistance(graph_t *g, uint_t i, uint_t j)
{
  double dist, xi, xj, yi, yj, zi, zj;
  slow_assert(!isArcOrEdge(g, i, j));

  xi  = g->contattr[g->x_index][i];
  xj  = g->contattr[g->x_index][j];
  yi  = g->contattr[g->y_index][i];
  yj  = g->contattr[g->y_index][j];
  zi  = g->contattr[g->z_index][i];
  zj  = g->contattr[g->z_index][j];
  

  if (isnan(xi) || isnan(yi) || isnan(zi) ||
      isnan(xj) || isnan(yj) || isnan(zj)) {
    return 0;
  }
  else {
    dist = euclidean_distance(xi, yi, zi, xj, yj, zj);
    return dist;
  }
}


/******************Attribute interaction (categorical) ***********************/

/* 
 * Change statistic for interaction effect of two categorical attributes 
 * matching: adding arc i->j increases statistic by 1 when
 * nodes i and j have the same categorical attribute value attribute a
 * and also have the same categorical attribute value for attribute b
 * (note a and b are different attributes, they don't have to have the same
 * value)
 */
double changeMatchingInteraction(graph_t *g, uint_t i, uint_t j,
                                 uint_t a, uint_t b)
{
  slow_assert(!isArcOrEdge(g, i, j));
  return g->catattr[a][i] != CAT_NA && g->catattr[a][j] != CAT_NA &&
    g->catattr[b][i] != CAT_NA && g->catattr[b][j] != CAT_NA &&
    g->catattr[a][i] == g->catattr[a][j] &&
    g->catattr[b][i] == g->catattr[b][j];
}



/*****************************************************************************
 *
 * other external functions
 *
 ****************************************************************************/

/*
 *
 * Compute the change statistics for addition of arc i->j
 * This involves summing over all the statistics specified for structural
 * effects, nodal attribute effects, dyadic covariate effects, and attribute
 * interaction effects.
 *
 * Parameters:
 *   g      - digraph object. Modifed if performMove is true.
 *   i      - node source of arc being added (or deleted)
 *   j      - node dest of arc being added (or deleted)
 *            (Note that the edge or arc i-j must NOT exist in g passed
 *            to this function - deletion is in the context of removing
 *            the edge before calling this [it will be added back in
 *            the sampler if the delete move is rejected])
 *   n      - number of parameters (length of theta vector and total
 *            number of change statistic functions)
 *   n_attr - number of attribute change stats functions
 *   n_dyadic -number of dyadic covariate change stats funcs
 *   n_attr_interaction - number of attribute interaction change stats funcs
 *   change_stats_funcs - array of pointers to change statistics functions
 *                        length is n-n_attr-n_dyadic-n_attr_interaction
 *   lambda_values      - array of lambda values for change stats funcs
 *   attr_change_stats_funcs - array of pointers to change statistics functions
 *                             length is n_attr
 *   dyadic_change_stats_funcs - array of pointers to dyadic change stats funcs
 *                             length is n_dyadic
 *   attr_interaction_change_stats_funcs - array of points to attribute 
 *                                        interaction change statstistics
 *                                        functinons. Length is
 *                                        n_attr_interaction.
 *   attr_indices   - array of n_attr attribute indices (index into g->binattr
 *                    or g->catattr) corresponding to attr_change_stats_funcs
 *                    E.g. for Sender effect on the first binary attribute,
 *                    attr_indices[x] = 0 and attr_change_stats_funcs[x] =
 *                    changeSender
 *   exponent_values    - array of exponent values for attr change stats funcs
 *                          length is n_attr 
 *   attr_interaction_pair_indices - array of n_attr_interaction attribute pair
 *                                   indices (as above, but each element is
 *                                   a pair of such indices) for attribute
 *                                   interaction effects.
 *   theta  - array of n parameter values corresponding to change stats funcs
 *   isDelete - TRUE if arc is being deleted (statistics negated then)
 *   changestats - (OUT) array of n change statistics values corresponding to
 *                 change stats funcs. Allocated by caller.
 *
 * Return value:
 *   Sum of all change statistics for additionof arc i->j
 */
double calcChangeStats(graph_t *g, uint_t i, uint_t j,
                       uint_t n, uint_t n_attr, uint_t n_dyadic,
                       uint_t n_attr_interaction,
                       change_stats_func_t *change_stats_funcs[],
                       double               lambda_values[],
                       attr_change_stats_func_t *attr_change_stats_funcs[],
                       dyadic_change_stats_func_t *dyadic_change_stats_funcs[],
                       attr_interaction_change_stats_func_t 
                                        *attr_interaction_change_stats_funcs[],
                       uint_t attr_indices[],
                       double exponent_values[],
                       uint_pair_t attr_interaction_pair_indices[],
                       const double theta[],
                       bool isDelete,
                       double changestats[])
{
  double total = 0;  /* sum of theta*changestats */
  uint_t l, param_i = 0;
  
  /* structural effects */
  for (l = 0; l < n - n_attr - n_dyadic - n_attr_interaction; l++) { 
    changestats[param_i] = (*change_stats_funcs[l])(g, i, j,
                                                    lambda_values[param_i]);
    total += theta[param_i] * (isDelete ? -1 : 1) * changestats[param_i];
    param_i++;
  }
  /* nodal attribute effects */
  for (l = 0; l < n_attr; l++) {
    changestats[param_i] = (*attr_change_stats_funcs[l])
      (g, i, j, attr_indices[l], isDelete, exponent_values[param_i]);
    total += theta[param_i] * (isDelete ? -1 : 1) * changestats[param_i];
    param_i++;
  }
  /* dyadic covariate effects */
  for (l = 0; l < n_dyadic; l++) {
    changestats[param_i] = (*dyadic_change_stats_funcs[l])(g, i, j);
    total += theta[param_i] * (isDelete ? -1 : 1) * changestats[param_i];
    param_i++;
  }
  /* attribute pair interaction effects */
  for (l = 0; l < n_attr_interaction; l++) {
    changestats[param_i] = (*attr_interaction_change_stats_funcs[l])
      (g, i, j, attr_interaction_pair_indices[l].first,
       attr_interaction_pair_indices[l].second);
    total += theta[param_i] * (isDelete ? -1 : 1) * changestats[param_i]; 
    param_i++;
  }
  return total;
}


/*
 *
 * Compute the observed statistics for the empty graph (no arcs; all nodes
 * are isolates). [The number of nodes is always fixed here].
 * 
 * This is zero for all of the statistics except Isolates, where it is
 * the number of nodes. (Any statistics added in the future where the value
 * for the empty graph is not zero will have to be added here).
 * 
 * Parameters:
 *   g      - digraph object.
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
 *   attr_interaction_change_stats_funcs - array of points to attribute 
 *                                        interaction change statstistics
 *                                        functinons. Length is
 *                                        n_attr_interaction.
 *   attr_indices   - array of n_attr attribute indices (index into g->binattr
 *                    or g->catattr) corresponding to attr_change_stats_funcs
 *                    E.g. for Sender effect on the first binary attribute,
 *                    attr_indices[x] = 0 and attr_change_stats_funcs[x] =
 *                    changeSender
 *   exponent_values    - array of exponent values for attr change stats funcs
 *                        same length as attr_change_stats_funcs
 *   attr_interaction_pair_indices - array of n_attr_interaction attribute pair
 *                                   indices (as above, but each element is
 *                                   a pair of such indices) for attribute
 *                                   interaction effects.
 *   emptystats - (OUT) array of n observed statistics values corresponding to
 *                 change stats funcs. Allocated by caller.
 *
 * Return value:
 *   Pointer to emptystats array (parameter)
 */
double *empty_graph_stats(graph_t *g,
                          uint_t n, uint_t n_attr, uint_t n_dyadic,
                          uint_t n_attr_interaction,
                          change_stats_func_t *change_stats_funcs[],
                          double lambda_values[],
                          attr_change_stats_func_t *attr_change_stats_funcs[],
                          dyadic_change_stats_func_t *dyadic_change_stats_funcs[],
                          attr_interaction_change_stats_func_t 
                          *attr_interaction_change_stats_funcs[],
                          uint_t attr_indices[],
                          double exponent_values[],                          
                          uint_pair_t attr_interaction_pair_indices[],
                          double emptystats[])
{
  uint_t l, param_i = 0;
  
  (void)attr_change_stats_funcs; /* unused parameter */
  (void)lambda_values;           /* unused parameter */
  (void)dyadic_change_stats_funcs; /* unused parameter */
  (void)exponent_values;         /* unused parameter */
  (void)attr_interaction_change_stats_funcs; /* unused parameter */
  (void)attr_indices; /* unused parameter */
  (void)attr_interaction_pair_indices; /* unused parameter */

  /* structural effects */
  for (l = 0; l < n - n_attr - n_dyadic - n_attr_interaction; l++) { 
    /* It is valid to compare function pointer this way in C */
    /* https://stackoverflow.com/questions/14985423/function-pointer-equality-in-c */
    if (change_stats_funcs[l] == changeIsolates) {
      /* The Isolates statistic is the number of nodes for empty graph */
      emptystats[param_i] = g->num_nodes;
    } else if (change_stats_funcs[l] == changeBipartiteIsolatesA) {
      emptystats[param_i] = g->num_A_nodes;
    } else if (change_stats_funcs[l] == changeBipartiteIsolatesB) {
      emptystats[param_i] = g->num_B_nodes;
    } else {
      emptystats[param_i] = 0;
    }
    param_i++;
  }
  /* nodal attribute effects */
  for (l = 0; l < n_attr; l++) {
    emptystats[param_i] = 0;
    param_i++;
  }
  /* dyadic covariate effects */
  for (l = 0; l < n_dyadic; l++) {
    emptystats[param_i] = 0;
    param_i++;
  }
  /* attribute pair interaction effects */
  for (l = 0; l < n_attr_interaction; l++) {
    emptystats[param_i] = 0;
    param_i++;
  }
  return emptystats;
}
