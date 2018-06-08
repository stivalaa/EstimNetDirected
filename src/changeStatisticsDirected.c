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
      delta += pow(1-1/lambda, g->mixTwoPathMatrix[INDEX2D(i, v, g->num_nodes)]);
  }
  for (k = 0; k < g->indegree[i]; k++) {
    v = g->revarclist[i][k];
    if (v == i || v == j)
      continue;
    if (isArc(g, v, j))
      delta += pow(1-1/lambda, g->mixTwoPathMatrix[INDEX2D(v, j, g->num_nodes)]);
  }
  delta += lambda * (1 - pow(1-1/lambda,
                             g->mixTwoPathMatrix[INDEX2D(i, j, g->num_nodes)]));
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
    assert(isArc(g, v, i));
    if (v == i || v == j)
      continue;
    if (isArc(g, j, v)) {
      delta +=
        pow(1-1/lambda, g->mixTwoPathMatrix[INDEX2D(i, v, g->num_nodes)]) +
        pow(1-1/lambda, g->mixTwoPathMatrix[INDEX2D(v, j, g->num_nodes)]);
    }
  }
  delta +=
    lambda * (1 - pow(1-1/lambda, 
                       g->mixTwoPathMatrix[INDEX2D(j, i, g->num_nodes)]));
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
    delta += pow(1-1/lambda, g->mixTwoPathMatrix[INDEX2D(i, v, g->num_nodes)]);
  }
  for (k = 0; k < g->indegree[i]; k++) {
    v = g->revarclist[i][k];
    if (v == i || v == j)
      continue;
    delta += pow(1-1/lambda, g->mixTwoPathMatrix[INDEX2D(v, j, g->num_nodes)]);
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
    delta += pow(1-1/lambda, g->outTwoPathMatrix[INDEX2D(j, v, g->num_nodes)]);
  }
  return delta;
}

/*
 * Change statisic for alternating two-paths A2P-TD (shared popularity +
 * muliptle two-paths), adjusting for multiple counting
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
