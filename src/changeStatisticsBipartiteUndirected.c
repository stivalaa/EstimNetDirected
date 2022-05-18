/*****************************************************************************
 * 
 * File:    changeStatisticsBipartiteUndirected.c
 * Author:  Alex Stivala
 * Created: May 2022
 *
 * Functions to compute graph change statistics for undirected
 * bipartite (two-mode) graphs. Each function takes a pointer to a
 * graph struct, and two node numbers i and j and returns the value of
 * the change statistic for adding the edge i -- j where i is a node
 * of MODE_A and j is a node of MODE_B.
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
 * Some of these functions are adapted from the original BPNet
 * code by Peng Wang:
 *
 *   Wang P, Robins G, Pattison P. PNet: A program for the simulation
 *   and estimation of exponential random graph models. Melbourne
 *   School of Psychological Science, The University of
 *   Melbourne. 2006. http://www.melnet.org.au/s/PNetManual.pdf
 *
 * And for the definitions of the change statistics:
 *
 *   Wang, P., Sharpe, K., Robins, G. L., & Pattison,
 *   P. E. (2009). Exponential random graph (p∗) models for affiliation
 *   networks. Social Networks, 31(1), 12-25.
 *
 *   Wang, P., Pattison, P., & Robins, G. (2013). Exponential random
 *   graph model specifications for bipartite networks—A dependence
 *   hierarchy. Social networks, 35(2), 211-222.
 * 
 * And also generally:
 * 
 *   Lusher, D., Koskinen, J., & Robins, G. (Eds.). (2013). Exponential
 *   random graph models for social networks: Theory, methods, and
 *   applications. New York, NY: Cambridge University Press.
 * 
 * especially Ch. 10:
 *
 *   Weng, P. (2013). Exponential random graph model extensions:
 *   models for multiple networks and bipartite networks. In
 *   "Exponential random graph models for social networks: Theory,
 *   methods, and applications." (pp. 115-129). New York, NY: Cambridge
 *   University Press.
 *
 * The reference for the MPNet software also defines parameters for
 * bipartite networks as a special case of multilevel networks:
 *
 *   Wang, P., Robins, G., Pattison, P., & Koskinen,
 *   J. (2014). Program for the Simulation and Estimation of (p*)
 *   Exponential Random Graph Models for Multilevel
 *   Networks. Melbourne School of Psychological Sciences, The
 *   University of Melbourne.
 *   https://www.melnet.org.au/s/MPNetManual.pdf
 *
 * As well as the statnet ergm terms, and references for specific
 * change statistics where indicated.
 *
 * Note that in the references above, the two modes are referred to as
 * 'A' and 'P' for 'association' and 'people' respectively, and
 * conventionally represented graphically (red) circle for 'P' nodes
 * and (blue) squares for 'A' nodes. This derives from the common use
 * of bipartite graphs to represent affiliation networks (e.g. company
 * directors on boards) where in the affiliation matrix, the people
 * (P) are rows and their affiliations (A) the columns. However here
 * we are going to use just name the node types A and B, where the
 * type A ndoes correspond to the 'A' in BPNet and type B to the 'P'
 * in BPNet [note this seems to be the other way around in the MPNet
 * naming convention according to the MPNet manual, as, confusingly,
 * in the MPNet manual, type 'A' nodes are shown as blue squares and
 * type 'B' as red circles.  Actually using MPNet shows that, by
 * comparing the observed statistics to confirm, XACB corresponds to
 * Kcp, XASA corresponds to KSa, and XASB corresponds to Ksp - so
 * the names used here are consistent with MPNet e.g. XACB in MPNet
 * is AltKCyclesB here].
 *
 * Do NOT compile with -ffast-math on gcc as we depend on IEEE handling of NaN
 *
 ****************************************************************************/

#include <math.h>
#include <assert.h>
#include "changeStatisticsBipartiteUndirected.h"


/*****************************************************************************
 *
 * local functions
 *
 ****************************************************************************/

/*
 * number of s-stars (s >=2) for a vertex v
 */
static ulong_t num_s_stars(const graph_t *g, uint_t v, ulong_t s)
{
  ulong_t d, num, i;
  ulong_t count = 0;
  assert(s >= 2);
  d = g->degree[v];
  if (s == 2) {
    return d * (d - 1) / 2;
  } else if (d >= s) {
    num = d;
    for (i = 1; i < s; i++) {
      num *= (d - i);
    }
    count += num / factorial(s);
    return count;
  }
  return count;
}

/*
 * change statistic for an s-star (s >= 2) for a vertex v
 */
static ulong_t change_s_stars(const graph_t *g, uint_t v, ulong_t s)
{
  assert(s >= 2);
  return s == 2 ? g->degree[v] : num_s_stars(g, v, s - 1);
}


/*****************************************************************************
 *
 * change statistics functions
 *
 ****************************************************************************/


/************************* Structural ****************************************/


/*
 * Change statistic for 2-stars for type A nodes
 */
double changeBipartiteTwoStarsA(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)lambda; /* unused parameters */
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  return change_s_stars(g, i, 2);
}

/*
 * Change statistic for 2-stars for type B nodes
 */
double changeBipartiteTwoStarsB(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)lambda; /* unused parameters */
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  return change_s_stars(g, j, 2);
}

/*
 * Change statistic for 3-stars for type A nodes
 */
double changeBipartiteThreeStarsA(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)lambda; /* unused parameters */
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  return change_s_stars(g, i, 3);
}

/*
 * Change statistic for 3-stars for type B nodes
 */
double changeBipartiteThreeStarsB(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)lambda; /* unused parameters */
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  return change_s_stars(g, j, 3);
}



/*
 * Change statistic for bipartite 4-cycle
 */
double changeBipartiteFourCycle(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t v,k;
  uint_t delta = 0;
  (void)lambda; /* unused parameters */  
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);

  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    if (v != i){
      delta += GET_A2PATH_ENTRY(g, v, i);
    }
  }
  return (double)delta;
}

/*
 * Change statistic for alternating k-stars for type A nodes
 */
double changeBipartiteAltStarsA(graph_t *g, uint_t i, uint_t j, double lambda)
{
  assert(lambda > 1);
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  return lambda * (1 - POW_LOOKUP(1-1/lambda, g->degree[i]));
}

/*
 * Change statistic for alternating k-stars for type B nodes
 */
double changeBipartiteAltStarsB(graph_t *g, uint_t i, uint_t j, double lambda)
{
  assert(lambda > 1);
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  return lambda * (1 - POW_LOOKUP(1-1/lambda, g->degree[j]));
}


/*
 * Change statistic for alternating k-cycles for type A nodes
 */
double changeBipartiteAltKCyclesA(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t k,v;
  double delta = 0;
  assert(lambda > 1);
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  for (k = 0; k < g->degree[i]; k++) {
    v = g->edgelist[i][k];
    if (v != j) {
      delta += POW_LOOKUP(1-1/lambda, GET_B2PATH_ENTRY(g, j, v));
    }
  }
  return delta;
}

/*
 * Change statistic for alternating k-cycles for type B nodes
 */
double changeBipartiteAltKCyclesB(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t k,v;
  double delta = 0;
  assert(lambda > 1);
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    if (v != i) {
      delta += POW_LOOKUP(1-1/lambda, GET_A2PATH_ENTRY(g, i, v));
    }
  }
  return delta;
}
