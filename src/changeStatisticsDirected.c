/*****************************************************************************
 * 
 * File:    changeStatisticsDirected.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * Functions to compute directed graph change statistics. Each
 * function takes a pointer to a digraph struct, and two node numbers
 * i and j and returns the value of the change statistic for adding
 * the arc i -> j (which must not already exist in the graph).
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
#include "changeStatisticsDirected.h"


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
 * Change statistic for Arc
 */
double changeArc(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)g; (void)i; (void)j; (void)lambda; /* unused parameters */
  assert(g->is_directed);
  return 1;
}

/*
 * Change statistic for Reciprocity
 */
double changeReciprocity(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)lambda; /* unused parameter */
  assert(g->is_directed);
  if (i == j) {
    return 0;
  } else {
    return isArc(g, j, i);
  }
}

/*
 * Change statistic for Sink 
 */
double changeSink(graph_t *g, uint_t i, uint_t j, double lambda)
{
  double delta = 0;
  (void)lambda; /* unused parameter */
  assert(g->is_directed);
  if (g->outdegree[i] == 0 && g->indegree[i] != 0) {
    delta--;
  }
  if (i != j && g->outdegree[j] == 0 && g->indegree[j] == 0) {
    delta++;
  }
  return delta;
}

/*
 * Change statistic for Source
 */
double changeSource(graph_t *g, uint_t i, uint_t j, double lambda)
{
  double delta = 0;
  (void)lambda; /* unused parameter */
  assert(g->is_directed);
  if (i != j && g->outdegree[i] == 0 && g->indegree[i] == 0) {
    delta++;
  }
  if (g->indegree[j] == 0 && g->outdegree[j] != 0) {
    delta--;
  }
  return delta;
}



/*
 * Change statistic for in-2-star (triad census 021U; but note that since
 * these staistics, unlike motifs, are not induced subgraphs, this also
 * counts, in some cases multiple times, 111D, 030T, 201, 120D, 120U,
 * 120C, 210, 300)
 */
double changeInTwoStars(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)lambda; /* unused parameter */
  (void)i; /* unused parameter */
  assert(g->is_directed);
  return g->indegree[j];
}

/*
 * Change statistic for out-2-star (triad census 021D; but not that since
 * these statistics, unlike motifs, are not induced subgraphs, this also counts,
 * in some cases multiple times,  111D, 111U, 030T, 201, 120D, 120U, 120C, 
 * 210, 300)
 */
double changeOutTwoStars(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)lambda; /* unused parameter */
  (void)j; /* unused parameter */
  return g->outdegree[i];
}

/*
 * Change statistic for transitive triangle (triad census 030T; but note
 * that since these statistics, unlike motifs, are not induced subgraphs,
 * this also coutns 120D, 120U, and 300).
 */
double changeTransitiveTriad(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t v,k,l,w;
  uint_t  delta = 0;
  (void)lambda; /* unused parameter */
  assert(g->is_directed);
  if (i == j) {
    return 0;
  }
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
double changeCyclicTriad(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t v,k;
  uint_t  delta = 0;
  (void)lambda; /* unused parameter */
  assert(g->is_directed);
  if (i == j) {
    return 0;
  }
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
double changeAltInStars(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t jindegree = g->indegree[j];
  (void)i; /*unused parameter*/
  assert(lambda > 1);
  assert(g->is_directed);
  return lambda * (1 - POW_LOOKUP(1-1/lambda, jindegree));
}

/*
 * Change statistic for alternating k-out-stars (activity spread, AoutS)
 */
double changeAltOutStars(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t ioutdegree = g->outdegree[i];
  (void)j;/*unused parameter*/
  assert(lambda > 1);
  assert(g->is_directed);
  return lambda * (1 - POW_LOOKUP(1-1/lambda, ioutdegree));
}

/*
 * Change statistic for alternating k-triangles AT-T (path closure)
 */
double changeAltKTrianglesT(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t v,k;
  double  delta = 0;
  assert(lambda > 1);
  assert(g->is_directed);
  
  if (i == j) {
    return 0;
  }

  for (k = 0; k < g->outdegree[i]; k++) {
    v = g->arclist[i][k];
    if (v == i || v == j)
      continue;
    if (isArc(g, j, v))
      delta += POW_LOOKUP(1-1/lambda,
                   GET_MIX2PATH_ENTRY(g, i, v));
  }
  for (k = 0; k < g->indegree[i]; k++) {
    v = g->revarclist[i][k];
    if (v == i || v == j)
      continue;
    if (isArc(g, v, j))
      delta += POW_LOOKUP(1-1/lambda,
                   GET_MIX2PATH_ENTRY(g, v, j));
  }
  delta += lambda * (1 - POW_LOOKUP(1-1/lambda,
                             GET_MIX2PATH_ENTRY(g, i, j)));
  return delta;
}


/*
 * Change statistic for alternating k-triangles AT-C (cyclic closure)
 */
double changeAltKTrianglesC(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t v,k;
  double delta =0;
  assert(lambda > 1);
  assert(g->is_directed);

  if (i == j) {
    return 0;
  }

  for (k = 0; k < g->indegree[i]; k++) {
    v = g->revarclist[i][k];
    slow_assert(isArc(g, v, i));
    if (v == i || v == j)
      continue;
    if (isArc(g, j, v)) {
      delta +=
        POW_LOOKUP(1-1/lambda, GET_MIX2PATH_ENTRY(g, i, v)) +
        POW_LOOKUP(1-1/lambda, GET_MIX2PATH_ENTRY(g, v, j));
    }
  }
  delta +=
    lambda * (1 - POW_LOOKUP(1-1/lambda, 
                      GET_MIX2PATH_ENTRY(g, j, i)));
  return delta;
}

/*
 * Change statistic for alternating k-triangles AT-D (popularity closure)
 */
double changeAltKTrianglesD(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t v,k;
  double delta = 0;
  assert(lambda > 1);
  assert(g->is_directed);

  if (i == j) {
    return 0;
  }

  for (k = 0; k < g->outdegree[i]; k++) {
    v = g->arclist[i][k];
    if (v == i || v == j)
      continue;
    if (isArc(g, j, v)) {
      delta +=
        POW_LOOKUP(1-1/lambda, GET_OUT2PATH_ENTRY(g, j, v));
    }
    if (isArc(g, v, j)) {
      delta += 
        POW_LOOKUP(1-1/lambda, GET_OUT2PATH_ENTRY(g, v, j));
    }
  }
  delta +=
    lambda * (1 - POW_LOOKUP(1-1/lambda, 
                      GET_OUT2PATH_ENTRY(g, i, j)));

  return delta;
}

/*
 * Change statistic for alternating k-triangles AT-U (activity closure)
 */
double changeAltKTrianglesU(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t v,k;
  double delta = 0;
  assert(lambda > 1);
  assert(g->is_directed);

  if (i == j) {
    return 0;
  }

  for (k = 0; k < g->indegree[j]; k++) {
    v = g->revarclist[j][k];
    if (v == i || v == j)
      continue;
    if (isArc(g, i, v)) {
      delta +=
        POW_LOOKUP(1-1/lambda, GET_IN2PATH_ENTRY(g, i, v));
    }
    if (isArc(g, v, i)) {
      delta += 
        POW_LOOKUP(1-1/lambda, GET_IN2PATH_ENTRY(g, v, i));
    }
  }
  delta +=
    lambda * (1 - POW_LOOKUP(1-1/lambda, 
                      GET_IN2PATH_ENTRY(g, i, j)));
  return delta;
}

/*
 * Change statistics for alternating two-path A2P-T (multiple 2-paths)
 */
double changeAltTwoPathsT(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t v,k;
  double delta = 0;
  assert(lambda > 1);
  assert(g->is_directed);

  if (i == j) {
    return 0;
  }

  for (k = 0; k < g->outdegree[j]; k++) {
    v = g->arclist[j][k];
    if (v == i || v == j)
      continue;
    delta += POW_LOOKUP(1-1/lambda, GET_MIX2PATH_ENTRY(g, i, v));
  }
  for (k = 0; k < g->indegree[i]; k++) {
    v = g->revarclist[i][k];
    if (v == i || v == j)
      continue;
    delta += POW_LOOKUP(1-1/lambda, GET_MIX2PATH_ENTRY(g, v, j));
  }

  return delta;
}

/*
 * Change statistic for alternating two-paths A2P-D (shared popularity) 
 */
double changeAltTwoPathsD(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t v,k;
  double delta = 0;
  assert(lambda > 1);
  assert(g->is_directed);

  if (i == j) {
    return 0;
  }

  for (k = 0; k < g->outdegree[i]; k++) {
    v = g->arclist[i][k];
    if (v == i || v == j) 
      continue;
    delta += POW_LOOKUP(1-1/lambda, GET_OUT2PATH_ENTRY(g, j, v));
  }
  return delta;
}

/*
 * Change statistic for alternating two-paths A2P-U (shared activity) 
 */
double changeAltTwoPathsU(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t v,k;
  double delta = 0;
  assert(lambda > 1);
  assert(g->is_directed);

  if (i == j) {
    return 0;
  }

  for (k = 0; k < g->indegree[j]; k++) {
    v = g->revarclist[j][k];
    if (v == i || v == j) 
      continue;
    delta += POW_LOOKUP(1-1/lambda, GET_IN2PATH_ENTRY(g, i, v));
  }

  return delta;
}

/*
 * Change statisic for alternating two-paths A2P-TD (shared popularity +
 * multiple two-paths), adjusting for multiple counting
 */
double changeAltTwoPathsTD(graph_t *g, uint_t i, uint_t j, double lambda)
{
  return 0.5 * (changeAltTwoPathsT(g, i, j, lambda) +
                changeAltTwoPathsD(g, i, j, lambda));
}


/*
 * Change statistic for Loop Interaction (number of arcs between
 * two nodes, both of which have a self-edge).
 * Note allowLoops = True must be set for this to work
 */
double changeLoopInteraction(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)lambda; /* unused parameter */
  uint_t delta = 0;
  uint_t k;
  assert(g->is_directed);

  /* first case, adding an arc between two nodes that both have self-edges */
  if (i != j && has_loop(g, i) && has_loop(g, j)) {
    ++delta;
  }
  else if (i == j) {
    /* second case, adding a self-edge to a node that has an arc to or from
       a node with a self-edge */
    for (k = 0; k < g->outdegree[i]; k++) {
      if (has_loop(g, g->arclist[i][k])) {
        ++delta;
      }
    }
    for (k = 0; k < g->indegree[i]; k++) {
      if (has_loop(g, g->revarclist[i][k])) {
        ++delta;
      }
    }
  }
  return (double)delta;
}

/************************* Actor attribute (binary) **************************/


/*
 * Change statistic for Sender
 */
double changeSender(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete) 
{
  (void)j; (void)isDelete; /*unused parameters*/
  assert(g->is_directed);
  return g->binattr[a][i] != BIN_NA && g->binattr[a][i];
}

/*
 * Change statistic for receiver
 */
double changeReceiver(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)i; (void)isDelete; /*unused parameters*/
  assert(g->is_directed);
  return g->binattr[a][j] != BIN_NA && g->binattr[a][j];
}


/********************* Actor attribute (categorical) *************************/


/*
 * Change statistic for categorical matching reciprocity
 */
double changeMatchingReciprocity(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)isDelete; /* unused parameter*/
  assert(g->is_directed);
  if (i == j) {
    return 0;
  } else {
    return g->catattr[a][i] != CAT_NA && g->catattr[a][j] != CAT_NA &&
           g->catattr[a][i] == g->catattr[a][j] && isArc(g, j, i);
  }
}


/*
 * Change statistic for categorical mismatching reciprocity
 */
double changeMismatchingReciprocity(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)isDelete; /* unused parameter */
  assert(g->is_directed);
  if (i == j) {
    return 0;
  } else {
    return g->catattr[a][i] != CAT_NA && g->catattr[a][j] != CAT_NA &&
           g->catattr[a][i] != g->catattr[a][j] && isArc(g, j, i);
  }
}


/*
 * Change statistic for different category transitive triangle:
 * transitive triangle i -> k, k -> j, i -> j where
 * i has a different value of the categorical attribute than that of j
 * and of k.
 *
 */
double changeMismatchingTransitiveTriad(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  uint_t v,k,l,w;
  uint_t  delta = 0;
  (void)isDelete; /* unused parameter */
  assert(g->is_directed);
  
  if (i == j) {
    return 0;
  }

  for (k = 0; k < g->outdegree[i]; k++) {
    v = g->arclist[i][k];
    if (v == i || v == j)
      continue;
    if (isArc(g, j, v) && g->catattr[a][i] != CAT_NA &&
        g->catattr[a][j] != CAT_NA && g->catattr[a][v] != CAT_NA &&
        g->catattr[a][i] != g->catattr[a][j] &&
        g->catattr[a][i] != g->catattr[a][v])
      delta++;
    if (isArc(g, v, j) && g->catattr[a][i] != CAT_NA &&
        g->catattr[a][j] != CAT_NA && g->catattr[a][v] != CAT_NA &&
        g->catattr[a][i] != g->catattr[a][j] &&
        g->catattr[a][i] != g->catattr[a][v])
      delta++;
  }
  for (l = 0; l < g->indegree[i]; l++) {
    w = g->revarclist[i][l];
    if (w == i || w == j)
      continue;
    if (isArc(g, w, j) && g->catattr[a][i] != CAT_NA &&
        g->catattr[a][j] != CAT_NA && g->catattr[a][w] != CAT_NA &&
        g->catattr[a][w] != g->catattr[a][i] &&
        g->catattr[a][w] != g->catattr[a][j])
      delta++;
  }
  return (double)delta;
}


/*
 * Change statistic for different category transitive ties
 * as defined in:
 *
 *   Schmid, C. S., Chen, T. H. Y., & Desmarais, B. A. (2021).
 *   Generative Dynamics of Supreme Court Citations:
 *   Analysis with a New Statistical Network Model. arXiv preprint
 *   arXiv:2101.07197
 *
 * where it is described as the "different term transitivity"
 * statistic ("difftransties" in cERGM), in the context of the
 * citation ERGM (cERGM) where the categorical attribute is the term
 * (time period).
 *
 * From the R manpage for difftransties in cERGM R package:
 * (https://github.com/desmarais-lab/Supreme_Court_Citation_Network)
 *
 *   'difftransties(attrname)' (directed) _Transitive ties with
 *        different sender attribute:_ This term adds one statistic,
 *        equal to the number of ties i-->j such that there exists a
 *        two-path from i to j, if for the nodal attribute 'attrname',
 *        the sending node i is different from j and the node on the
 *        two-path.
 *
 * (Note that this is a variation on the "transitiveties" statistic in
 * statnet, which allows a version where all three nodes involved must
 * have the same value of an attribute to be counted; here instead
 * the node that has two outgoing ties must have a different value from
 * the two nodes to which it sends those ties.)
 *
 * See d_difftransties in changestats.scc.c from cERGM (and c_transitivities
 * in changestats.c in ergm).
 *
 * This change statistic function is (so far) the only one whih uses
 * the isDelete parameter: instead of always computing the change statistic
 * for the arc i->j being added (the return value being negated by
 * the caller, calcChangeStats()), this function will temporarily
 * insert the arc i->j if isDelete is true, so that code similar to
 * the original implementation d_difftransties in cERGM can be used).
 * Note that unlike statnet change statistic functions however we don't
 * negate the value for deleting a tie here, since it is always done
 * in calcChangeStats().
 */
double changeMismatchingTransitiveTies(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  uint_t u, v,k,l;
  int L2th, L2tu, L2uh;
  int delta = 0;
  int ochange = isDelete ? -1 : 0;
  assert(g->is_directed);
  
  if (i == j) {
    return 0;
  }

  if (isDelete) {
    insertArc(g, i, j); /* temporarily put arc i->j in for delete */
  }

  L2th = 0;
  for (k = 0; k < g->outdegree[j]; k++) {
    u = g->arclist[j][k];
    if (isArc(g, i, u) && g->catattr[a][i] != CAT_NA &&
        g->catattr[a][j] != CAT_NA && g->catattr[a][u] != CAT_NA &&
        g->catattr[a][i] != g->catattr[a][j] &&
        g->catattr[a][i] != g->catattr[a][u]) {
      L2tu = ochange;
      for (l = 0; l < g->indegree[u]; l++) {
        v = g->revarclist[u][l];
        if (isArc(g, i, v) && g->catattr[a][i] != g->catattr[a][v]) {
          L2tu++;
          if (L2tu >0)
            break;
        }
      }
      delta += (L2tu == 0);
    }
  }
  for (k = 0; k < g->indegree[j]; k++) {
    u = g->revarclist[j][k];
    if (isArc(g, i, u) && g->catattr[a][i] != CAT_NA &&
        g->catattr[a][j] != CAT_NA && g->catattr[a][u] != CAT_NA &&
        g->catattr[a][i] != g->catattr[a][j] &&
        g->catattr[a][i] != g->catattr[a][u]) {
      L2th++;
    }
    if (isArc(g, u, i) && g->catattr[a][i] != CAT_NA &&
        g->catattr[a][j] != CAT_NA && g->catattr[a][u] != CAT_NA &&
        g->catattr[a][i] != g->catattr[a][u] &&
        g->catattr[a][j] != g->catattr[a][u]) {
      L2uh = ochange;
      for (l = 0; l < g->outdegree[u]; l++) {
        v = g->arclist[u][l];
        if (isArc(g, v, j) && g->catattr[a][v] != g->catattr[a][u]) {
          L2uh++;
          if (L2uh > 0)
            break;
        }
      }
      delta += (L2uh == 0);
    }
  }
  delta += (L2th > 0);

  if (isDelete) {
    removeArc(g, i, j); /* remove temporary i->j arc added at start */
  }

  return (double)delta;
}

/********************* Actor attribute (continuous) *************************/


/*
 * Change statistic for continuous Sender
 */
double changeContinuousSender(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)j; (void)isDelete; /*unused parameters*/
  assert(g->is_directed);
  if (isnan(g->contattr[a][i]))
    return 0;
  else
    return g->contattr[a][i];
}

/*
 * Change statistic for continuous Receiver
 */
double changeContinuousReceiver(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)i; (void)isDelete; /*unused parameters*/
  assert(g->is_directed);
  if (isnan(g->contattr[a][j]))
    return 0;
  else
    return g->contattr[a][j];
}



/*
 * Change statistic for continuous diff (absolute difference of attribute)
 * reciprocity
 */
double changeDiffReciprocity(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)isDelete; /* unused parameter */
  assert(g->is_directed);
  if (i == j)
    return 0;
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
double changeDiffSign(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)isDelete; /* unused parameter */
  assert(g->is_directed);
  if (isnan(g->contattr[a][i]) || isnan(g->contattr[a][j]))
    return 0;
  else
    return signum(g->contattr[a][i] - g->contattr[a][j]);
}


/*
 * Change statistic for signed continuous difference
 * for attr_i - attr_j if attr_i > attr_j and zero otherwise.
 * (so larger as sending node has higher attribute value and 0
 * when receiving node has higher attribute value, and the value is always
 * positive).
 * Like diff(dir="t-h", sign.action="posonly") in statnet.
 * ("tail" is sender, "head" is receiver)
 */
double changeDiffDirSR(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)isDelete; /* unused parameter */
  assert(g->is_directed);
  if (isnan(g->contattr[a][i]) || isnan(g->contattr[a][j]))
  {
    return 0;
  }
  else {
    if (g->contattr[a][i] > g->contattr[a][j])
      return g->contattr[a][i] - g->contattr[a][j];
    else
      return 0;
  }
}

/*
 * Change statistic for signed continuous difference
 * for attr_j - attr_i if attr_j > attr_i and zero otherwise.
 * (so larger as recieving node has higher attribute value and 0
 * when sending node has higher attribute value, and the value is always
 * positive).
 * Like diff(dir="h-t", sign.action="posonly") in statnet.
 * ("tail" is sender, "head" is receiver)
 */
double changeDiffDirRS(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)isDelete; /* unused parameter */
  assert(g->is_directed);
  if (isnan(g->contattr[a][i]) || isnan(g->contattr[a][j]))
  {
    return 0;
  }
  else {
    if (g->contattr[a][j] > g->contattr[a][i])
      return g->contattr[a][j] - g->contattr[a][i];
    else
      return 0;
  }
}





