/*****************************************************************************
 * 
 * File:    changeStatisticsBipartiteUndirected.c
 * Author:  Alex Stivala
 * Created: May 2022
 *
 * Functions to compute graph change statistics for undirected
 * bipartite (two-mode) graphs. Each function takes a pointer to a
 * graph struct, and two node numbers i and j and returns the value of
 * the change statistic for adding the edge i -- j (which must not
 * already exist in the graph) where i is a node of MODE_A and j is a
 * node of MODE_B.
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
 *   Wang, P. (2013). Exponential random graph model extensions:
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
 * b1nodematch and b2nodematch (statnet ergm names) are defined in:
 *
 *    Bomiriya, R. P. (2014). Topics in exponential random graph
 *    modeling. (Doctoral dissertation, Pennsylvania State University).
 *    https://etda.libraries.psu.edu/files/final_submissions/9857
 *
 *    Bomiriya, R. P., Kuvelkar, A. R., Hunter, D. R., & Triebel,
 *    S. (2023). Modeling Homophily in Exponential-Family Random Graph
 *    Models for Bipartite Networks. arXiv preprint
 *    arXiv:2312.05673. https://arxiv.org/abs/2312.05673
 *
 * Do NOT compile with -ffast-math on gcc as we depend on IEEE handling of NaN
 *
 ****************************************************************************/

#include <math.h>
#include <assert.h>
#include "changeStatisticsUndirected.h"
#include "changeStatisticsBipartiteUndirected.h"


/*****************************************************************************
 *
 * utility functions
 *
 ****************************************************************************/

/*
 * Compute value of integer x raised to the power y (double in [0,1]),
 * defining pow(0, 0) = 0 as per Bomiryia et al. (2023) [see p. 7
 * after eqn (7)]
 */
static double pow0(uint_t x, double y)
{
  return (x == 0 && DOUBLE_APPROX_EQ(y, 0)) ? 0 : pow(x, y);
}

/*
 * In Bomiriya et al. (2023) [Appendix A] change statistics definition
 * for node-centered (alpha-based homompily), t(i, j, k) in eqn (12)
 * is defined as the number of two-paths from i to j not passing
 * through k. Since only one two-path from i to j can pass through k,
 * this is just the number of two-paths from i to j, less one if
 * k is a neighbour of both i and j.
 */
static uint_t twopaths_not_via_k_A(const graph_t *g, uint_t i, uint_t j, uint_t k)
{
  uint_t count = GET_A2PATH_ENTRY(g, i, j);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_A);
  assert(bipartite_node_mode(g, k) == MODE_B);
  if (isEdge(g, i, k) && isEdge(g, j, k)) {
    count--;
  }
  return count;
}

static uint_t twopaths_not_via_k_B(const graph_t *g, uint_t i, uint_t j, uint_t k)
{
  uint_t count = GET_B2PATH_ENTRY(g, i, j);
  assert(bipartite_node_mode(g, i) == MODE_B);
  assert(bipartite_node_mode(g, j) == MODE_B);
  assert(bipartite_node_mode(g, k) == MODE_A);
  if (isEdge(g, i, k) && isEdge(g, j, k)) {
    count--;
  }
  return count;
}


/*****************************************************************************
 *
 * change statistics functions
 *
 ****************************************************************************/


/************************* Structural ****************************************/


/*
 * Change statistic for 2-stars for type A nodes
 * (Sa2 in BPNet, XStar2A in MPNet)
 */
double changeBipartiteTwoStarsA(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)lambda; /* unused parameters */
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  return (double)change_s_stars(g, i, 2);
}

/*
 * Change statistic for 2-stars for type B nodes
 * (Sp2 in BPNet, XStar2B in MPNet)
 */
double changeBipartiteTwoStarsB(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)lambda; /* unused parameters */
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  return (double)change_s_stars(g, j, 2);
}

/*
 * Change statistic for 3-stars for type A nodes
 * (Sa3 in BPNet, XStar3A in MPNet)
 */
double changeBipartiteThreeStarsA(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)lambda; /* unused parameters */
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  return (double)change_s_stars(g, i, 3);
}

/*
 * Change statistic for 3-stars for type B nodes
 * (Sp3 in BPNet, XStar3B in MPNet)
 */
double changeBipartiteThreeStarsB(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)lambda; /* unused parameters */
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  return (double)change_s_stars(g, j, 3);
}


/*
 * Change statistic for alternating k-stars for type A nodes
 * (K-Sa in BPNet, ASA in MPNet)
 */
double changeBipartiteAltStarsA(graph_t *g, uint_t i, uint_t j, double lambda)
{
  assert(lambda > 1);
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  return lambda * (1 - POW_LOOKUP(1-1/lambda, g->degree[i]));
}

/*
 * Change statistic for alternating k-stars for type B nodes
 * (K-Sp in BPNet, ASB in MPNet)
 */
double changeBipartiteAltStarsB(graph_t *g, uint_t i, uint_t j, double lambda)
{
  assert(lambda > 1);
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  return lambda * (1 - POW_LOOKUP(1-1/lambda, g->degree[j]));
}


/*
 * Change statistic for alternating k-cycles for type A nodes
 * (K-Ca in BPNet, XACA in MPNet)
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
  slow_assert(!isEdge(g, i, j));
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
 * (K-Cp in BPNet, XACB in MPNet)
 * 
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
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    if (v != i) {
      delta += POW_LOOKUP(1-1/lambda, GET_A2PATH_ENTRY(g, i, v));
    }
  }
  return delta;
}



/*
 * Change statistic for alternating k-4-cycles for type A nodes
 * (new change statistic suggested in email (basically paper outline, 
 * with spreadsheet attachments for literature search, examples, etc.)
 * "idea for a (slightly) new bipartite change statistic" sent 23 Nov 2022)
 */
double changeBipartiteAltK4CyclesA(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t k,v;
  double delta = 0;
  assert(lambda > 1);
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[i]; k++) {
    v = g->edgelist[i][k];
    if (v != j) {
      uint_t b2p = GET_B2PATH_ENTRY(g, j, v);
      delta += POW_LOOKUP(1-1/lambda, MAX(b2p - 1, 0));
    }
  }
  return delta;
}

/*
 * Change statistic for alternating k-4-cycles for type B nodes
 * (new change statistic suggested in email (basically paper outline, 
 * with spreadsheet attachments for literature search, examples, etc.)
 * "idea for a (slightly) new bipartite change statistic" sent 23 Nov 2022)
 */
double changeBipartiteAltK4CyclesB(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t k,v;
  double delta = 0;
  assert(lambda > 1);
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    if (v != i) {
      uint_t a2p = GET_A2PATH_ENTRY(g, i, v);
      delta += POW_LOOKUP(1-1/lambda, MAX(a2p - 1, 0));
    }
  }
  return delta;
}

/*
 * Change statistic for Isolates for mode A nodes
 */
double changeBipartiteIsolatesA(graph_t *g, uint_t i, uint_t j, double lambda)
{
  double delta = 0;
  (void)lambda; /* unused parameter */
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  if (g->degree[i] == 0) {
    delta--;
  }
  return delta;
}

/*
 * Change statistic for Isolates for mode B nodes
 */
double changeBipartiteIsolatesB(graph_t *g, uint_t i, uint_t j, double lambda)
{
  double delta = 0;
  (void)lambda; /* unused parameter */
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  if (g->degree[j] == 0) {
    delta--;
  }
  return delta;
}

/************************* Actor attribute (binary) **************************/


/*
 * Change statistic for Bipartite Activity for type A nodes
 * (RA in BPNet, XEdgeA (binary) in MPNet)
 */
double changeBipartiteActivityA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)isDelete; /*unused parameters*/
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  return (g->binattr[a][i] == BIN_NA ? 0 : g->binattr[a][i]);
}

/*
 * Change statistic for Bipartite Activity for type B nodes
 * (RP in BPNet, XEdgeB (binary) in MPNet)
 */
double changeBipartiteActivityB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)isDelete; /*unused parameters*/
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  return (g->binattr[a][j] == BIN_NA ? 0 : g->binattr[a][j]);
}


/*
 * Change statistic for Bipartite Interaction
 * (rAP in BPNet, ??? in MPNet) - not needed, just use changeInteraction()
 */


/*********************** Actor attribute (continuous) ************************/

/*
 * Change statistic for Bipartite Continuous Activity for type A nodes
 * (RAC in BPNet, XEdgeA (continuous) in MPNet)
 */
double changeBipartiteContinuousActivityA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)isDelete; /*unused parameters*/
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  if (isnan(g->contattr[a][i])) {
    return 0;
  } else {
    return g->contattr[a][i];
  }
}

/*
 * Change statistic for Bipartite Continuous Activity for type B nodes
 * (RPC in BPNet, XEdgeB (continuous) in MPNet)
 */
double changeBipartiteContinuousActivityB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)isDelete; /*unused parameters*/
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  if (isnan(g->contattr[a][j])) {
    return 0;
  } else {
    return g->contattr[a][j];
  }
}


/*
 * Change statistic for Bipartite 2-path sum for type A nodes
 * (TSOACS in BPNet, X2StarASum in MPNet)
 */
double changeBipartiteTwoPathSumA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  uint_t k,v;
  double delta = 0;
  (void)isDelete; /*unused parameters*/
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    assert(v != j);
    if (v != i) {
      if (!isnan(g->contattr[a][i]) &&
          !isnan(g->contattr[a][v])) {
        delta += g->contattr[a][i] + g->contattr[a][v];
      }
    }
  }
  return (double)delta;
}


/*
 * Change statistic for Bipartite 2-path sum for type B nodes
 * (TSOPCS in BPNet, X2StarBSum in MPNet)
 */
double changeBipartiteTwoPathSumB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  uint_t k,v;
  double delta = 0;
  (void)isDelete; /*unused parameters*/
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[i]; k++) {
    v = g->edgelist[i][k];
    assert(v != i);
    if (v != j) {
      if (!isnan(g->contattr[a][j]) &&
          !isnan(g->contattr[a][v])) {
        delta += g->contattr[a][j] + g->contattr[a][v];
      }
    }
  }
  return (double)delta;
}

/*
 * Change statistic for Bipartite 2-path abs difference for type A nodes
 * (TSOACD in BPNet, X2StarADifference in MPNet)
 */
double changeBipartiteTwoPathDiffA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  uint_t k,v;
  double delta = 0;
  (void)isDelete; /*unused parameters*/
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    assert(v != j);
    if (v != i) {
      if (!isnan(g->contattr[a][i]) &&
          !isnan(g->contattr[a][v])) {
        delta += fabs(g->contattr[a][i] - g->contattr[a][v]);
      }
    }
  }
  return (double)delta;
}


/*
 * Change statistic for Bipartite 2-path abs diff for type B nodes
 * (TSOPCD in BPNet, X2StarBDifference in MPNet)
 */
double changeBipartiteTwoPathDiffB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  uint_t k,v;
  double delta = 0;
  (void)isDelete; /*unused parameters*/
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[i]; k++) {
    v = g->edgelist[i][k];
    assert(v != i);
    if (v != j) {
      if (!isnan(g->contattr[a][j]) &&
          !isnan(g->contattr[a][v])) {
        delta += fabs(g->contattr[a][j] - g->contattr[a][v]);
      }
    }
  }
  return (double)delta;
}


/*
 * Change statistic for Bipartite Continuous Sum
 * (RAPC in BPNet, XEdgeABSum in MPNet) - not needed, just use changeSum()
 */

/*********************** Actor attribute (categorical) ************************/

/*
 * Change statistic for Bipartite 2-path matching for type A nodes
 * (2path_match_A in BPNet, X2StarAMatch in MPNet)
 */
double changeBipartiteTwoPathMatchingA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  uint_t k,v;
  uint_t delta = 0;
  (void)isDelete; /*unused parameters*/
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    assert(v != j);
    if (v != i) {
      if (g->catattr[a][i] != CAT_NA &&
          g->catattr[a][v] != CAT_NA &&
          g->catattr[a][i] == g->catattr[a][v]) {
        delta++;
      }
    }
  }
  return (double)delta;
}


/*
 * Change statistic for Bipartite 2-path mismatching for type B nodes
 * (2path_mismatch_P in BPNet, X2StarBMismatch in MPNet)
 */
double changeBipartiteTwoPathMatchingB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  uint_t k,v;
  uint_t delta = 0;
  (void)isDelete; /*unused parameters*/
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[i]; k++) {
    v = g->edgelist[i][k];
    assert(v != i);
    if (v != j) {
      if (g->catattr[a][j] != CAT_NA &&
          g->catattr[a][v] != CAT_NA &&
          g->catattr[a][j] == g->catattr[a][v]) {
        delta++;
      }
    }
  }
  return (double)delta;
}

/*
 * Change statistic for Bipartite 2-path mismatching for type A nodes
 * (2path_mismatch_A in BPNet, X2StarAMismatch in MPNet)
 */
double changeBipartiteTwoPathMismatchingA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  uint_t k,v;
  uint_t delta = 0;
  (void)isDelete; /*unused parameters*/
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    assert(v != j);
    if (v != i) {
      if (g->catattr[a][i] != CAT_NA &&
          g->catattr[a][v] != CAT_NA &&
          g->catattr[a][i] != g->catattr[a][v]) {
        delta++;
      }
    }
  }
  return (double)delta;
}


/*
 * Change statistic for Bipartite 2-path matching for type B nodes
 * (2path_match_P in BPNet, X2StarBMatch in MPNet)
 */
double changeBipartiteTwoPathMismatchingB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  uint_t k,v;
  uint_t delta = 0;
  (void)isDelete; /*unused parameters*/
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[i]; k++) {
    v = g->edgelist[i][k];
    assert(v != i);
    if (v != j) {
      if (g->catattr[a][j] != CAT_NA &&
          g->catattr[a][v] != CAT_NA &&
          g->catattr[a][j] != g->catattr[a][v]) {
        delta++;
      }
    }
  }
  return (double)delta;
}


/**************** Actor attribute (categorical, exponent) ********************/

/*
 * Change statistic for Bipartite node-centered (alpha-based) homophily
 * for type A node (b1nodematch(alpha) statnet ergm term)
 *
 * alpha is the exponent in the range [0, 1]
 *
 * b1nodematch and b2nodematch (statnet ergm names) are defined in:
 *
 *  Bomiriya, R. P. (2014). Topics in exponential random graph
 *  modeling. (Doctoral dissertation, Pennsylvania State University).
 *  https://etda.libraries.psu.edu/files/final_submissions/9857
 *
 *  Bomiriya, R. P., Kuvelkar, A. R., Hunter, D. R., & Triebel,
 *  S. (2023). Modeling Homophily in Exponential-Family Random Graph
 *  Models for Bipartite Networks. arXiv preprint
 *  arXiv:2312.05673. https://arxiv.org/abs/2312.05673
 *
 * This change statistic is defined by equation (12) in Bomiriya et al. (2023)
 *
 */
double changeBipartiteNodematchAlphaA(graph_t *g, uint_t i, uint_t j, uint_t a, double alpha)
{
  uint_t k, v;
  double delta = 0;

  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[i]; k++) {
    v = g->edgelist[i][k];
    assert(v != i);
    if (v != j) {
      if (g->catattr[a][j] != CAT_NA &&
          g->catattr[a][v] != CAT_NA &&
          g->catattr[a][j] != g->catattr[a][v]) {
        /* Note pow0 defines pow0(0, 0) = 0
           as per Bomiryia et al. (2023) [see p. 7 after eqn (7)] */
        delta += (pow0(twopaths_not_via_k_A(g, i, v, j) + 1, alpha) -
                  pow0(twopaths_not_via_k_A(g, i, v, j), alpha));
      }
    }
  }
  return delta;
}

/*
 * Change statistic for Bipartite edge-centered (beta-based) homophily
 * for type A node (b1nodematch(beta) statnet ergm term)
 *
 * beta is the exponent in the range [0, 1]
 *
 * b1nodematch and b2nodematch (statnet ergm names) are defined in:
 *
 *  Bomiriya, R. P. (2014). Topics in exponential random graph
 *  modeling. (Doctoral dissertation, Pennsylvania State University).
 *  https://etda.libraries.psu.edu/files/final_submissions/9857
 *
 *  Bomiriya, R. P., Kuvelkar, A. R., Hunter, D. R., & Triebel,
 *  S. (2023). Modeling Homophily in Exponential-Family Random Graph
 *  Models for Bipartite Networks. arXiv preprint
 *  arXiv:2312.05673. https://arxiv.org/abs/2312.05673
 *
 * This change statistic is defined by equation (14) in Bomiriya et al. (2023)
 */
double changeBipartiteNodematchBetaA(graph_t *g, uint_t i, uint_t j, uint_t a, double beta)
{
  uint_t k, v;
  uint_t u = 0; /* number of edges to j from nodes (not i) matching i */
  double delta = 0;

  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));

  /* count u = number of edges to j from nodes (not i) matching i */
  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    assert(v != j);
    if (v != i) {
      if (g->catattr[a][i] != CAT_NA &&
          g->catattr[a][v] != CAT_NA &&
          g->catattr[a][i] != g->catattr[a][v]) {
        u++;
      }
    }
  }
  /* Note pow0 defines pow0(0, 0) = 0 as per Bomiryia et al. (2023)
     [see p. 7 after eqn (7)] */
  delta = 0.5 * ( (1 + u)*pow0(u, beta) - u*pow0(u - 1, beta) );
  return delta;
}
