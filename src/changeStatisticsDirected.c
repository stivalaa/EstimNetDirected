/*****************************************************************************
 * 
 * File:    changeStatisticsDirected.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * Functions to compute directed graph change statistics. Each
 * function takes a pointer to a digraph struct, and two node numbers
 * i and j and returns the value of the change statistic for adding
 * the arc i -> j.
 *
 * For change statistics dependent on a nodal attribute, there is
 * an additional parameter a which is the index of the attribute
 * to use.
 *
 * These functions are adapted from the original PNet code by Peng Wang:
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
 ****************************************************************************/

#include <math.h>
#include <assert.h>
#include "changeStatisticsDirected.h"

   
/*****************************************************************************
 *
 * constants
 *
 ****************************************************************************/

/* 
 * lambda decay parameter for alternating statistics 
 */
const double lambda = 2.0; /* TODO make it a configuration setting */

/*****************************************************************************
 *
 * local functions
 *
 ****************************************************************************/

/*
 * signum function, returns -1 for negative x, +1 for positive x, else 0
 */
static double signum(double x)
{
  /* https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c */
  return (0 < x) - (x < 0);
}



/*****************************************************************************
 *
 * change statistics functions
 *
 ****************************************************************************/


/************************* Structural ****************************************/


/* 
 * Change statistic for Ac
 */
double changeArc(const digraph_t *g, uint_t i, uint_t j)
{
  (void)g; (void)i; (void)j; /* unused parameters */
  return 1;
}

/*
 * Change statistic for Reciprocity
 */
double changeReciprocity(const digraph_t *g, uint_t i, uint_t j)
{
  return isArc(g, j, i);
}

/*
 * Change statistic for Sink 
 */
double changeSink(const digraph_t *g, uint_t i, uint_t j)
{
  double delta = 0;
  if (g->outdegree[i] == 0 && g->indegree[i] != 0) {
    delta--;
  }
  if (g->outdegree[j] == 0 && g->indegree[j] == 0) {
    delta++;
  }
  return delta;
}

/*
 * Change statistic for Source
 */
double changeSource(const digraph_t *g, uint_t i, uint_t j)
{
  double delta = 0;
  if (g->outdegree[i] == 0 && g->indegree[i] == 0) {
    delta++;
  }
  if (g->indegree[j] == 0 && g->outdegree[j] != 0) {
    delta--;
  }
  return delta;
}

/*
 * Change statistic for Isolates
 */
double changeIsolates(const digraph_t *g, uint_t i, uint_t j)
{
  double delta = 0;
  if (g->indegree[i] == 0 && g->outdegree[i] == 0) {
    delta--;
  }
  if (g->indegree[j] == 0 && g->outdegree[j] == 0) {
    delta--;
  }
  return delta;
}

/*
 * Change statistic for two-path (triad census 021C)
 * also known as TwoMixStar
 */
double changeTwoPath(const digraph_t *g, uint_t i, uint_t j)
{
  return g->indegree[i] + g->outdegree[j] - (isArc(g, j, i) ? 2 : 0);
}

/*
 * Change statistic for in-2-star (triad census 021U)
 */
double changeInTwoStars(const digraph_t *g, uint_t i, uint_t j)
{
  (void)i; /* unused parameter */
  return g->indegree[j];
}

/*
 * Change statistic for out-2-star (triad census 021D)
 */
double changeOutTwoStars(const digraph_t *g, uint_t i, uint_t j)
{
  (void)j; /* unused parameter */
  return g->outdegree[i];
}

/*
 * Change statistic for transitive triangle (triad census 030T)
 */
double changeTransitiveTriad(const digraph_t *g, uint_t i, uint_t j)
{
  uint_t v,k,l,w;
  uint_t  delta = 0;
  for (k = 0; k < g->outdegree[i]; k++) {
    v = g->arclist[i][k];
    if (v == i || v == j)
      continue;
    if (isArc(g, j, v))
      delta++;
    if (isArc(g, v, j))
      delta++;
  }
  for (l = 0; l < g->indegree[i]; l++) {
    w = g->revarclist[i][l];
    if (w == i || w == j)
      continue;
    if (isArc(g, w, j))
      delta++;
  }
  return (double)delta;
}

/*
 * Change statistic for cyclic triangle (triad census 030C)
 */
double changeCyclicTriad(const digraph_t *g, uint_t i, uint_t j)
{
  uint_t v,k;
  uint_t  delta = 0;
  for (k = 0; k < g->indegree[i]; k++) {
    v = g->revarclist[i][k];
    if (v == i || v == j)
      continue;
    if (isArc(g, j, v))
      delta++;
  }
  return (double)delta;
}

/*
 * Change statistic for alternating k-in-stars (popularity spread, AinS)
 */
double changeAltInStars(const digraph_t *g, uint_t i, uint_t j)
{
  uint_t jindegree = g->indegree[j];
  (void)i; /*unused parameter*/
  assert(lambda > 1);
  return lambda * (1 - pow(1-1/lambda, jindegree));
}

/*
 * Change statistic for alternating k-out-stars (activity spread, AoutS)
 */
double changeAltOutStars(const digraph_t *g, uint_t i, uint_t j)
{
  uint_t ioutdegree = g->outdegree[i];
  (void)j;/*unused parameter*/
  assert(lambda > 1);
  return lambda * (1 - pow(1-1/lambda, ioutdegree));
}

/*
 * Change statistic for alternating k-triangles AT-T (path closure)
 */
double changeAltKTrianglesT(const digraph_t *g, uint_t i, uint_t j)
{
  uint_t v,k;
  double  delta = 0;
  assert(lambda > 1);
  for (k = 0; k < g->outdegree[i]; k++) {
    v = g->arclist[i][k];
    if (v == i || v == j)
      continue;
    if (isArc(g, j, v))
      delta += pow(1-1/lambda,
                   GET_MIX2PATH_ENTRY(g, i, v));
  }
  for (k = 0; k < g->indegree[i]; k++) {
    v = g->revarclist[i][k];
    if (v == i || v == j)
      continue;
    if (isArc(g, v, j))
      delta += pow(1-1/lambda,
                   GET_MIX2PATH_ENTRY(g, v, j));
  }
  delta += lambda * (1 - pow(1-1/lambda,
                             GET_MIX2PATH_ENTRY(g, i, j)));
  return delta;
}

/*
 * Change statistic for alternating k-triangles AT-C (cyclic closure)
 */
double changeAltKTrianglesC(const digraph_t *g, uint_t i, uint_t j)
{
  uint_t v,k;
  double delta =0;
  assert(lambda > 1);
  for (k = 0; k < g->indegree[i]; k++) {
    v = g->revarclist[i][k];
    /*removed as slows significantly: assert(isArc(g, v, i));*/
    if (v == i || v == j)
      continue;
    if (isArc(g, j, v)) {
      delta +=
        pow(1-1/lambda, GET_MIX2PATH_ENTRY(g, i, v)) +
        pow(1-1/lambda, GET_MIX2PATH_ENTRY(g, v, j));
    }
  }
  delta +=
    lambda * (1 - pow(1-1/lambda, 
                      GET_MIX2PATH_ENTRY(g, j, i)));
  return delta;
}

/*
 * Change statistic for alternating k-triangles AT-D (popularity closure)
 */
double changeAltKTrianglesD(const digraph_t *g, uint_t i, uint_t j)
{
  uint_t v,k;
  double delta = 0;
  assert(lambda > 1);
  for (k = 0; k < g->outdegree[i]; k++) {
    v = g->arclist[i][k];
    if (v == i || v == j)
      continue;
    if (isArc(g, j, v)) {
      delta +=
        pow(1-1/lambda, GET_OUT2PATH_ENTRY(g, j, v));
    }
    if (isArc(g, v, j)) {
      delta += 
        pow(1-1/lambda, GET_OUT2PATH_ENTRY(g, v, j));
    }
  }
  delta +=
    lambda * (1 - pow(1-1/lambda, 
                      GET_OUT2PATH_ENTRY(g, i, j)));
  return delta;
}

/*
 * Change statistic for alternating k-triangles AT-U (activity closure)
 */
double changeAltKTrianglesU(const digraph_t *g, uint_t i, uint_t j)
{
  uint_t v,k;
  double delta = 0;
  assert(lambda > 1);
  for (k = 0; k < g->indegree[j]; k++) {
    v = g->revarclist[j][k];
    if (v == i || v == j)
      continue;
    if (isArc(g, i, v)) {
      delta +=
        pow(1-1/lambda, GET_IN2PATH_ENTRY(g, i, v));
    }
    if (isArc(g, v, i)) {
      delta += 
        pow(1-1/lambda, GET_IN2PATH_ENTRY(g, v, i));
    }
  }
  delta +=
    lambda * (1 - pow(1-1/lambda, 
                      GET_IN2PATH_ENTRY(g, i, j)));
  return delta;
}

/*
 * Change statistics for alternating two-path A2P-T (multiple 2-paths)
 */
double changeAltTwoPathsT(const digraph_t *g, uint_t i, uint_t j)
{
  uint_t v,k;
  double delta = 0;
  assert(lambda > 1);
  for (k = 0; k < g->outdegree[j]; k++) {
    v = g->arclist[j][k];
    if (v == i || v == j)
      continue;
    delta += pow(1-1/lambda, GET_MIX2PATH_ENTRY(g, i, v));
  }
  for (k = 0; k < g->indegree[i]; k++) {
    v = g->revarclist[i][k];
    if (v == i || v == j)
      continue;
    delta += pow(1-1/lambda, GET_MIX2PATH_ENTRY(g, v, j));
  }
  return delta;
}

/*
 * Change statistic for alternating two-paths A2P-D (shared popularity) 
 */
double changeAltTwoPathsD(const digraph_t *g, uint_t i, uint_t j)
{
  uint_t v,k;
  double delta = 0;
  assert(lambda > 1);
  for (k = 0; k < g->outdegree[i]; k++) {
    v = g->arclist[i][k];
    if (v == i || v == j) 
      continue;
    delta += pow(1-1/lambda, GET_OUT2PATH_ENTRY(g, j, v));
  }
  return delta;
}

/*
 * Change statistic for alternating two-paths A2P-U (shared activity) 
 */
double changeAltTwoPathsU(const digraph_t *g, uint_t i, uint_t j)
{
  uint_t v,k;
  double delta = 0;
  assert(lambda > 1);
  for (k = 0; k < g->indegree[j]; k++) {
    v = g->revarclist[j][k];
    if (v == i || v == j) 
      continue;
    delta += pow(1-1/lambda, GET_IN2PATH_ENTRY(g, i, v));
  }
  return delta;
}

/*
 * Change statisic for alternating two-paths A2P-TD (shared popularity +
 * multiple two-paths), adjusting for multiple counting
 */
double changeAltTwoPathsTD(const digraph_t *g, uint_t i, uint_t j)
{
  return 0.5 * (changeAltTwoPathsT(g, i, j) + changeAltTwoPathsD(g, i, j));
}


/************************* Actor attribute (binary) **************************/


/*
 * Change statistic for Sender
 */
double changeSender(const digraph_t *g, uint_t i, uint_t j, uint_t a) 
{
  (void)j;/*unused parameter*/
  return g->binattr[a][i] != BIN_NA && g->binattr[a][i];
}

/*
 * Change statistic for receiver
 */
double changeReceiver(const digraph_t *g, uint_t i, uint_t j, uint_t a)
{
  (void)i;/*unused parameter*/
  return g->binattr[a][j] != BIN_NA && g->binattr[a][j];
}

/*
 * Change statistic for Interaction
 */
double changeInteraction(const digraph_t *g, uint_t i, uint_t j, uint_t a)
{
  return g->binattr[a][i] != BIN_NA && g->binattr[a][j] != BIN_NA &&
         g->binattr[a][i] && g->binattr[a][j];
}

/********************* Actor attribute (categorical) *************************/


/*
 * Change statistic for categorical matching
 */
double changeMatching(const digraph_t *g, uint_t i, uint_t j, uint_t a)
{
  return g->catattr[a][i] != CAT_NA && g->catattr[a][j] != CAT_NA &&
         g->catattr[a][i] == g->catattr[a][j];
}

/*
 * Change statistic for categorical matching reciprocity
 */
double changeMatchingReciprocity(const digraph_t *g, uint_t i, uint_t j, uint_t a)
{
  return g->catattr[a][i] != CAT_NA && g->catattr[a][j] != CAT_NA &&
         g->catattr[a][i] == g->catattr[a][j] && isArc(g, j, i);
}

/*
 * Change statistic for categorical mismatching
 */
double changeMismatching(const digraph_t *g, uint_t i, uint_t j, uint_t a)
{
  return g->catattr[a][i] != CAT_NA && g->catattr[a][j] != CAT_NA &&
         g->catattr[a][i] != g->catattr[a][j];
}

/*
 * Change statistic for categorical mismatching reciprocity
 */
double changeMismatchingReciprocity(const digraph_t *g, uint_t i, uint_t j, uint_t a)
{
  return g->catattr[a][i] != CAT_NA && g->catattr[a][j] != CAT_NA &&
         g->catattr[a][i] != g->catattr[a][j] && isArc(g, j, i);
}

/********************* Actor attribute (continuous) *************************/


/*
 * Change statistic for continuous Sender
 */
double changeContinuousSender(const digraph_t *g, uint_t i, uint_t j, uint_t a)
{
  (void)j;/*unused parameter*/
  if (isnan(g->contattr[a][i]))
    return 0;
  else
    return g->contattr[a][i];
}

/*
 * Change statistic for continuous Receiver
 */
double changeContinuousReceiver(const digraph_t *g, uint_t i, uint_t j, uint_t a)
{
  (void)i;/*unused parameter*/
  if (isnan(g->contattr[a][j]))
    return 0;
  else
    return g->contattr[a][j];
}


/*
 * Change statistic for continuous diff (absolute difference of attribute)
 */
double changeDiff(const digraph_t *g, uint_t i, uint_t j, uint_t a)
{
  if (isnan(g->contattr[a][i]) || isnan(g->contattr[a][j]))
    return 0;
  else
    return fabs(g->contattr[a][i] - g->contattr[a][j]);
}


/*
 * Change statistic for continuous diff (absolute difference of attribute)
 * reciprocity
 */
double changeDiffReciprocity(const digraph_t *g, uint_t i, uint_t j, uint_t a)
{
  if (isnan(g->contattr[a][i]) || isnan(g->contattr[a][j]))
    return 0;
  else  
    return fabs(g->contattr[a][i] - g->contattr[a][j]) * isArc(g, j, i);
}


/*
 * Change statistic for continuous diff sign (sign of difference of attribute)
 * for attr_i - attr_j (so +1 when sending node has higher attribute value and -1
 * when receiving node has higher attribute value).
 */
double changeDiffSign(const digraph_t *g, uint_t i, uint_t j, uint_t a)
{
  if (isnan(g->contattr[a][i]) || isnan(g->contattr[a][j]))
    return 0;
  else
    return signum(g->contattr[a][i] - g->contattr[a][j]);
}


/********************* Dyadic covariate (continuous) *************************/

/* 
 * Change steatistic for geographical distance between two nodes,
 * using for each node the pair of continuous attributes labelled 
 * as being latitude and longitude
 */
double changeGeoDistance(const digraph_t *g, uint_t i, uint_t j)
{
  double dist, lati, longi, latj, longj;

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
 * Change steatistic for logarithm of geographical distance between two nodes,
 * using for each node the pair of continuous attributes labelled 
 * as being latitude and longitude
 */
double changeLogGeoDistance(const digraph_t *g, uint_t i, uint_t j)
{
  double dist;

  dist = changeGeoDistance(g, i, j);
  if (dist > 0) {
    return log(dist);
  } else {
    return 0;
  }
}


/******************Attribute interaction (categorical) ***********************/

double changeMatchingInteraction(const digraph_t *g, uint_t i, uint_t j,
                                 uint_t a, uint_t b)
{
  /* TODO*/ return 0;
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
 *   n      - number of parameters (length of theta vector and total
 *            number of change statistic functions)
 *   n_attr - number of attribute change stats functions
 *   n_dyadic -number of dyadic covariate change stats funcs
 *   change_stats_funcs - array of pointers to change statistics functions
 *                        length is n-n_attr
 *   attr_change_stats_funcs - array of pointers to change statistics functions
 *                             length is n_attr
 *   dyadic_change_stats_funcs - array of pointers to dyadic change stats funcs
 *                             length is n_dyadic
 *   attr_indices   - array of n_attr attribute indices (index into g->binattr
 *                    or g->catattr) corresponding to attr_change_stats_funcs
 *                    E.g. for Sender effect on the first binary attribute,
 *                    attr_indices[x] = 0 and attr_change_stats_funcs[x] =
 *                    changeSender
 *   theta  - array of n parameter values corresponding to change stats funcs
 *   isDelete - TRUE if arc is being deleted (statistics negated then)
 *   changestats - (OUT) array of n change statistics values corresponding to
 *                 change stats funcs. Allocated by caller.
 *
 * Return value:
 *   Sum of all change statistics for additionof arc i->j
 */
double calcChangeStats(const digraph_t *g, uint_t i, uint_t j,
                       uint_t n, uint_t n_attr, uint_t n_dyadic,
                       change_stats_func_t *change_stats_funcs[],
                       attr_change_stats_func_t *attr_change_stats_funcs[],
                       dyadic_change_stats_func_t *dyadic_change_stats_funcs[],
                       uint_t attr_indices[], const double theta[],
                       bool isDelete,
                       double changestats[])
{
  double total = 0;  /* sum of theta*changestats */
  uint_t l, param_i = 0;
  
  /* structural effects */
  for (l = 0; l < n - n_attr - n_dyadic; l++) { 
    changestats[param_i] = (*change_stats_funcs[l])(g, i, j);
    total += theta[param_i] * (isDelete ? -1 : 1) * changestats[param_i];
    param_i++;
  }
  /* nodal attribute effects */
  for (l = 0; l < n_attr; l++) {
    changestats[param_i] = (*attr_change_stats_funcs[l])
      (g, i, j, attr_indices[l]);
    total += theta[param_i] * (isDelete ? -1 : 1) * changestats[param_i];
    param_i++;
  }
  /* dyadic covariate effects */
  for (l = 0; l < n_dyadic; l++) {
    changestats[param_i] = (*dyadic_change_stats_funcs[l])(g, i, j);
    total += theta[param_i] * (isDelete ? -1 : 1) * changestats[param_i];
    param_i++;
  }
  return total;
}
