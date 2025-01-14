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
 * type A nodss correspond to the 'A' in BPNet and type B to the 'P'
 * in BPNet [note this seems to be the other way around in the MPNet
 * naming convention according to the MPNet manual, as, confusingly,
 * in the MPNet manual, type 'A' nodes are shown as blue squares and
 * type 'B' as red circles.  Actually using MPNet shows that, by
 * comparing the observed statistics to confirm, XACB corresponds to
 * Kcp, XASA corresponds to KSa, and XASB corresponds to Ksp - so
 * the names used here are consistent with MPNet e.g. XACB in MPNet
 * is AltKCyclesB here].
 *
 * This is very confusing, so perhaps an example is necessary.
 * For this purpose I will use the Inouye-Pyke pollinator web data
 * from ../examples/bipartite/inouye_pyke_pollinators/.
 * This network has 91 pollinators and 42 plants.
 *
 * For BPNet, supplying an affiliation matrix with 91 rows and 42 columns,
 * we have from BPNet (in file start-statistics-inouye_pyke.txt):
 *
 *  ****This graph contains:****
 * vertices_A	91
 * vertices_P	42
 * L	281
 * Ksa	305.594483
 * Ksp	415.502266
 * Kcp	1167.875000
 *
 * (note this is already confusing, as although the supplied matrix
 * has 91 rows representing the pollinators ("Persons") and 42 columns
 * representing plants ("Affiliations") as conventional, BPNet takes the
 * rows as type A and the columns as type P).
 *
 * MPNet using the same input matrix gives (with 91 A nodes and
 * 42 B nodes in GUI under "Number of nodes", since A must be rows
 * and B columns), at the top of the setting file inouye_mpnet.pnet:
 *
 * Session_Name
 * inouye_mpnet
 * Session_Folder
 * C:\Users\user\switchdrive\Institution\USI\shared\ERGMXL\example_bipartite_networks\MPNet_estimations\ERGM\Inouye_Pyke_pollinator_web
 * Num_nodes_A
 * 91
 * Num_nodes_P
 * 42
 *
 * indicating that A in BPNet corresponds to A in MPNet and
 * P in BPNet corresponds to B in MPNet. And with this
 * model (inouye_mpnet_est.txt):
 *
 * Effects	Lambda	Parameter	Stderr	t-ratio	SACF
 * XEdge	2.0000	-5.5726	0.197	0.024	0.678	*
 * XASA	2.0000	0.7998	0.127	0.022	0.623	*
 * XASB	2.0000	0.7301	0.174	0.020	0.649	*
 * XACB	2.0000	0.0730	0.002	-0.012	0.806	*
 *
 * the observed statistics are (inouye_mpnet_est.txt):
 *
 * Observed graph statistics:
 * 281.00	305.59	415.50	1167.88	
 *
 * i.e. XEdge = 281, XASA = 305.39, XASB = 415.50, XACB = 1167.88,
 * showing that KSa in BPNet is XASA in MPNet, KSp in BPNet is XASB in
 * MPNET, and Kcp in BPNet is XACB in MPNet. This is also confirmed
 * with the GoF output (inouye_mpnet_gof.txt):
 *
 * Statistics      Observed        Mean    StdDev  t-ratio
 * XEdge   281.00000000    177.69300000    85.37181473     1.21008321
 * X3Path  15919.00000000  3548.88900000   3965.69826496   3.11927690      #
 * X4Cycle 830.00000000    84.73700000     111.69013310    6.67259479      #
 * XASA    305.59448338    163.21189795    112.98077897    1.26023724
 * XASB    415.50226605    231.82008951    146.54563247    1.25341283
 * XACA    586.65625000    224.56791016    163.67465093    2.21224446      #
 * XACB    1167.87500000   516.95787500    411.77137307    1.58077314
 *
 *
 * With EstimNetDirected (see ../examples/bipartite/inouye_pyke_pollinators/)
 * we have (stdout):
 *
 * Two-mode Graph with 133 vertices (91 mode A, 42 mode B) and 281 edges (density 0.0735217) [loops not allowed]
 *
 * and (obs_stats_ifd_inouye_pyke_pollinators_altkcycles_0.txt):
 *
 * Edge BipartiteAltStarsA(2) BipartiteAltStarsB(2) BipartiteAltKCyclesA(2) BipartiteAltKCyclesB(2)
 * 281 305.594 415.502 586.656 1167.88
 *
 *
 * so BipartiteAltStarsA   = XASA = Ksa
 *    BipartiteAltStasrB   = XASB = Ksp
 *    BipartiteAltKCyclesA = XACA = Kca
 *    BipartiteAltKCyclesB = XACB = Kcp
 *
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
double pow0(uint_t x, double y)
{
  return (x == 0 && DOUBLE_APPROX_EQ(y, 0)) ? 0 : pow(x, y);
}

/*****************************************************************************
 *
 * local utility functions
 *
 ****************************************************************************/


/*
 * In Bomiriya et al. (2023) [Appendix A] change statistics definition
 * for node-centered (alpha-based homompily), t(i, j, k) in eqn (12)
 * is defined as the number of two-paths from i to j not passing
 * through k. Since only one two-path from i to j can pass through k,
 * this is just the number of two-paths from i to j, less one if
 * k is a neighbour of both i and j.
 */
static uint_t twopaths_not_via_k(const graph_t *g, uint_t i, uint_t j,
				 uint_t k, bipartite_node_mode_e mode)
{
  uint_t count;
  if (mode == MODE_A) {
    count = GET_A2PATH_ENTRY(g, i, j);
  } else {
    count = GET_B2PATH_ENTRY(g, i, j);
  }
  assert(bipartite_node_mode(g, i) == mode);
  assert(bipartite_node_mode(g, j) == mode);
  assert(bipartite_node_mode(g, k) == other_mode(mode));
  if (isEdge(g, i, k) && isEdge(g, j, k)) {
    assert(count > 0); /* because there is a path from i to j via k */
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
 * (K-Ca in BPNet, XACA in MPNet) defined by eqn (6.14) in:
 *
 *   Wang, P., Sharpe, K., Robins, G. L., & Pattison,
 *   P. E. (2009). Exponential random graph (p∗) models for affiliation
 *   networks. Social Networks, 31(1), 12-25.
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
    assert(bipartite_node_mode(g, v) == MODE_B);
    if (v != j) {
      delta += POW_LOOKUP(1-1/lambda, GET_B2PATH_ENTRY(g, j, v));
    }
  }
  return delta;
}

/*
 * Change statistic for alternating k-cycles for type B nodes
 * (K-Cp in BPNet, XACB in MPNet) defined by eqn (6.14) in:
 *
 *   Wang, P., Sharpe, K., Robins, G. L., & Pattison,
 *   P. E. (2009). Exponential random graph (p∗) models for affiliation
 *   networks. Social Networks, 31(1), 12-25.
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
    assert(bipartite_node_mode(g, v) == MODE_A);
    if (v != i) {
      delta += POW_LOOKUP(1-1/lambda, GET_A2PATH_ENTRY(g, i, v));
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
double changeBipartiteActivityA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
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
double changeBipartiteActivityB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
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


/*
 * Change statistic for bipartite exactly one neighbour with binary attribute a
 * for type A nodes.
 *
 * The statistic counts the number of type A nodes that have exactly one
 * neighbour (therefore of type B) with the binary attribute a.
 *
 * Note that binary attribute a here is a binary attribute for type B nodes.
 */
double changeBipartiteExactlyOneNeighbourA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  uint_t num_neighbours_with_a = 0;
  uint_t delta = 0;
  uint_t k, v;
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[i]; k++) {
    v = g->edgelist[i][k];
    if (g->binattr[a][v] != BIN_NA && g->binattr[a][v]) {
      num_neighbours_with_a++;
      /* Note could shortcut and break out of loop as soon as
         num_neighbours_with_a == 2 as only need to know if 0, 1, or > 1
         but why complicate things? */
    }
  }
  /* the statistic can only change if j has binary attribute a */
  if (g->binattr[a][j] != BIN_NA && g->binattr[a][j]) {
    if (num_neighbours_with_a == 0) {
      /* if i has no neighbours with a and j has a, then i--j creates
       * a type A node with exactly one neihbour with a */
      delta++;
    } else if (num_neighbours_with_a == 1) {
      /* if i has exactly one neighbour with a, and j has a, then i--j
       * decreases by one the number of type A nodes with exactly one
       * neighbour with a */
      delta--;
    }
    /* if i has > 1 neighgours with a, no change in statistic */
  }
  return (double)delta;
}

/*
 * Change statistic for bipartite exactly one neighbour with binary attribute a
 * for type B nodes.
 *
 * The statistic counts the number of type B nodes that have exactly one
 * neighbour (therefore of type A) with the binary attribute a.
 *
 * Note that binary attribute a here is a binary attribute for type A nodes.
 */
double changeBipartiteExactlyOneNeighbourB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  uint_t num_neighbours_with_a = 0;
  uint_t delta = 0;
  uint_t k, v;
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    if (g->binattr[a][v] != BIN_NA && g->binattr[a][v]) {
      num_neighbours_with_a++;
      /* Note could shortcut and break out of loop as soon as
         num_neighbours_with_a == 2 as only need to know if 0, 1, or > 1
         but why complicate things? */
    }
  }
  /* the statistic can only change if i has binary attribute a */
  if (g->binattr[a][i] != BIN_NA && g->binattr[a][i]) {
    if (num_neighbours_with_a == 0) {
      /* if j has no neighbours with a and i has a, then i--j creates
       * a type B node with exactly one neihbour with a */
      delta++;
    } else if (num_neighbours_with_a == 1) {
      /* if j has exactly one neighbour with a, and i has a, then i--j
       * decreases by one the number of type B nodes with exactly one
       * neighbour with a */
      delta--;
    }
    /* if i has > 1 neighgours with a, no change in statistic */
  }
  return (double)delta;
}


/*********************** Actor attribute (continuous) ************************/

/*
 * Change statistic for Bipartite Continuous Activity for type A nodes
 * (RAC in BPNet, XEdgeA (continuous) in MPNet)
 */
double changeBipartiteContinuousActivityA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
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
double changeBipartiteContinuousActivityB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
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
double changeBipartiteTwoPathSumA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  uint_t k,v;
  double delta = 0;
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
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
double changeBipartiteTwoPathSumB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  uint_t k,v;
  double delta = 0;
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
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
double changeBipartiteTwoPathDiffA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  uint_t k,v;
  double delta = 0;
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
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
double changeBipartiteTwoPathDiffB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  uint_t k,v;
  double delta = 0;
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
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




/*********************** Actor attribute (categorical) ***********************/

/*
 * Change statistic for Bipartite 2-path matching for type A nodes
 * (2path_match_A in BPNet, X2StarAMatch in MPNet)
 */
double changeBipartiteTwoPathMatchingA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  uint_t k,v;
  uint_t delta = 0;
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
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
 * Change statistic for Bipartite 2-path matching for type B nodes
 * (2path_match_P in BPNet, X2StarBMatch in MPNet)
 */
double changeBipartiteTwoPathMatchingB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  uint_t k,v;
  uint_t delta = 0;
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
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
double changeBipartiteTwoPathMismatchingA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  uint_t k,v;
  uint_t delta = 0;
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
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
 * Change statistic for Bipartite 2-path mismatching for type B nodes
 * (2path_mismatch_P in BPNet, X2StarBMismatch in MPNet)
 */
double changeBipartiteTwoPathMismatchingB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  uint_t k,v;
  uint_t delta = 0;
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
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
 * for type A or B node (b1nodematch(alpha) or b2nodematch(alpha)
 * statnet ergm term). An extra parameter, mode, is passed which determines
 * if it is for A or B type node - the actual change statistic function with
 * the usual signature calls this with appropriate mode.
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
static double changeBipartiteNodematchAlpha(graph_t *g, uint_t i, uint_t j, uint_t a, double alpha, bipartite_node_mode_e mode)
{
  uint_t k, v;
  double delta = 0;
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == mode);
  assert(bipartite_node_mode(g, j) == other_mode(mode));
  slow_assert(!isEdge(g, i, j));
  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    assert(v != j);
    assert(bipartite_node_mode(g, v) == mode);
    if (v != i) {
      if (g->catattr[a][i] != CAT_NA &&
          g->catattr[a][v] != CAT_NA &&
          g->catattr[a][i] == g->catattr[a][v]) {
        /* Note pow0 defines pow0(0, 0) = 0
           as per Bomiryia et al. (2023) [see p. 7 after eqn (7)] */
        delta += (pow0(twopaths_not_via_k(g, i, v, j, mode) + 1, alpha) -
                  pow0(twopaths_not_via_k(g, i, v, j, mode), alpha));
      }
    }
  }
  return delta;
}

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
double changeBipartiteNodematchAlphaA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double alpha)
{
  (void)isDelete; /*unused parameter*/
  return changeBipartiteNodematchAlpha(g, i, j, a, alpha, MODE_A);
}

/*
 * Change statistic for Bipartite node-centered (alpha-based) homophily
 * for type B node (b2nodematch(alpha) statnet ergm term)
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
double changeBipartiteNodematchAlphaB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double alpha)
{
  (void)isDelete; /*unused parameter*/
  return changeBipartiteNodematchAlpha(g, j, i, a, alpha, MODE_B);
}

/*
 * Change statistic for Bipartite edge-centered (beta-based) homophily
 * for type A or B node (b1nodematch(beta) or b2nodemach(beta) statnet
 * ergm term). An extra parameter, mode, is passed which determines
 * if it is for A or B type node - the actual change statistic
 * function with the usual signature calls this with appropriate mode.
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
static double changeBipartiteNodematchBeta(graph_t *g, uint_t i, uint_t j, uint_t a, double beta, bipartite_node_mode_e mode)
{
  uint_t k, v;
  uint_t u = 0; /* number of edges to j from nodes (not i) matching i */
  double delta = 0;
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == mode);
  assert(bipartite_node_mode(g, j) == other_mode(mode));
  slow_assert(!isEdge(g, i, j));

  /* count u = number of edges to j from nodes (not i) matching i */
  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    assert(v != j);
    assert(bipartite_node_mode(g, v) == mode);
    if (v != i) {
      if (g->catattr[a][i] != CAT_NA &&
          g->catattr[a][v] != CAT_NA &&
          g->catattr[a][i] == g->catattr[a][v]) {
        u++;
      }
    }
  }
  /* Note pow0 defines pow0(0, 0) = 0 as per Bomiryia et al. (2023)
     [see p. 7 after eqn (7)] */
  delta = 0.5 * ( (1 + u)*pow0(u, beta) - u*pow0(u - 1, beta) );
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
double changeBipartiteNodematchBetaA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double beta)
{
  (void)isDelete; /*unused parameter*/
  return changeBipartiteNodematchBeta(g, i, j, a, beta, MODE_A);
}

/*
 * Change statistic for Bipartite edge-centered (beta-based) homophily
 * for type B node (b2nodematch(beta) statnet ergm term)
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
double changeBipartiteNodematchBetaB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double beta)
{
  (void)isDelete; /*unused parameter*/
  return changeBipartiteNodematchBeta(g, j, i, a, beta, MODE_B);
}


/*****************************************************************************
 *
 * experimental change statistics functions
 *
 ****************************************************************************/


/*
 * Change statistic for alternating k-4-cycles for type A nodes
 * (new change statistic suggested in email (basically paper outline, 
 * with spreadsheet attachments for literature search, examples, etc.)
 * "idea for a (slightly) new bipartite change statistic" sent 23 Nov 2022)
 *
 */
double changeBipartiteAltK4CyclesA(graph_t *g, uint_t i, uint_t j, double lambda)
{
  assert(lambda > 1);
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  return -1.0*(changeBipartiteAltKCyclesA(g, i, j, lambda) - change_s_stars(g, i, 2));
}

/*
 * Change statistic for alternating k-4-cycles for type B nodes
 * (new change statistic suggested in email (basically paper outline, 
 * with spreadsheet attachments for literature search, examples, etc.)
 * "idea for a (slightly) new bipartite change statistic" sent 23 Nov 2022)
 *
 */
double changeBipartiteAltK4CyclesB(graph_t *g, uint_t i, uint_t j, double lambda)
{
  assert(lambda > 1);
  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == MODE_A);
  assert(bipartite_node_mode(g, j) == MODE_B);
  slow_assert(!isEdge(g, i, j));
  return -1.0*(changeBipartiteAltKCyclesB(g, i, j, lambda) - change_s_stars(g, j, 2));
}



/*
 * Change statistic for number of 4-cycles at each node raised to a
 * power. The lambda parameter (> 1.0) (mis)used to specify the value
 * 1/lambda as the epxonent. Note this is not the same meaning of
 * lambda as its original use in the "alternating" parameters.
 *
 * An extra parameter, mode, is passed which determines
 * if it is for A or B type node - the actual change statistic
 * function with the usual signature calls this with appropriate mode.
 *
 */
static double changeBipartitePowerFourCycles(graph_t *g, uint_t i, uint_t j, double lambda, bipartite_node_mode_e mode)
{
  uint_t  v,k;
  ulong_t delta = 0;
  ulong_t vcount = 0;
  double  alpha = 1/lambda;
  double  change = 0;

  assert(g->is_bipartite);
  assert(!g->is_directed);
  assert(bipartite_node_mode(g, i) == mode);
  assert(bipartite_node_mode(g, j) == other_mode(mode));
  slow_assert(!isEdge(g, i, j));

  /* Number of four-cycles that node of selected mode is already inolved in */
  ulong_t count =  (bipartite_node_mode(g, i) == mode ?
                    num_four_cycles_node(g, i) :
                    num_four_cycles_node(g, j));
    
  /* change statistic for four-cycles */
  delta = changeFourCycles(g, i, j, lambda);
  change = pow(count + delta, alpha) - pow(count, alpha);

  /* neighbours of node of the opposite of the selected mode, so these
     are nodes of the selected mode */
  uint_t oppnode   = bipartite_node_mode(g, i) == mode ? j : i;
  uint_t othernode = bipartite_node_mode(g, i) == mode ? i : j;
  for (k = 0; k < g->degree[oppnode]; k++) {
    v = g->edgelist[oppnode][k];
    assert(bipartite_node_mode(g, v) == mode);
    vcount = num_four_cycles_node(g, v);
    if (bipartite_node_mode(g, v) == MODE_A) {
      delta = GET_A2PATH_ENTRY(g, v, othernode);
    } else {
      delta = GET_B2PATH_ENTRY(g, v, othernode);      
    }
    change += pow(vcount + delta, alpha) - pow(vcount, alpha);
  }
  return change;
}

/*
 * Change statistic for number of 4-cycles at each node raised to a
 * power. The lambda parameter (> 1.0) (mis)used to specify the value
 * 1/lambda as the epxonent. Note this is not the same meaning of
 * lambda as its original use in the "alternating" parameters.
 *
 * This change statistic for type A nodes.
 *
 */
double changeBipartitePowerFourCyclesA(graph_t *g, uint_t i, uint_t j, double lambda)
{
  return changeBipartitePowerFourCycles(g, i, j, lambda, MODE_A);
}

/*
 * Change statistic for number of 4-cycles at each node raised to a
 * power. The lambda parameter (> 1.0) (mis)used to specify the value
 * 1/lambda as the epxonent. Note this is not the same meaning of
 * lambda as its original use in the "alternating" parameters.
 *
 * This change statistic for type B nodes.
 *
 */
double changeBipartitePowerFourCyclesB(graph_t *g, uint_t i, uint_t j, double lambda)
{
  return changeBipartitePowerFourCycles(g, j, i, lambda, MODE_B);
}
