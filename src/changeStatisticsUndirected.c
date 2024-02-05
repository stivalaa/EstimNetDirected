/*****************************************************************************
 * 
 * File:    changeStatisticsUndirected.c
 * Author:  Alex Stivala
 * Created: January 2022
 *
 * Functions to compute graph change statistics for undirected graphs. Each
 * function takes a pointer to a graph struct, and two node numbers
 * i and j and returns the value of the change statistic for adding
 * the edge i -- j (which must not already exist in the graph).
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
#include "changeStatisticsUndirected.h"


/*****************************************************************************
 *
 * utility functions
 *
 ****************************************************************************/

/*
 * number of s-stars (s >=2) for a vertex v
 */
ulong_t num_s_stars(const graph_t *g, uint_t v, ulong_t s)
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
ulong_t change_s_stars(const graph_t *g, uint_t v, ulong_t s)
{
  assert(s >= 2);
  return s == 2 ? g->degree[v] : num_s_stars(g, v, s - 1);
}



/*
 * Binomial coefficient n choose 2
 */
static ulong_t n_choose_2(uint_t n)
{
  if (n < 2) {
    return 0;
  }
  return n * (n - 1) / 2;
}

/*
 * Count number of four-cycles that a particular node u is involved in.
 *
 * Note can also be used as for bipartite networks.
 */
uint_t num_four_cycles_node(const graph_t *g, uint_t u)
{
  uint_t k,v,i,j,l;
  uint_t count = 0;
  bool *visited = safe_calloc(g->num_nodes, sizeof(bool));

  if (g->is_bipartite) {
    /* iterate over all nodes that are distance 2 from u */
    for (k = 0; k < g->degree[u]; k++){
      i = g->edgelist[u][k];
      assert(bipartite_node_mode(g, i) != bipartite_node_mode(g, u));
      for (l = 0; l < g->degree[i]; l++) {
        j = g->edgelist[i][l];
        assert(bipartite_node_mode(g, j) == bipartite_node_mode(g, u));
        if (j != u && !visited[j]) {
          visited[j] = TRUE;
          if (bipartite_node_mode(g, u) == MODE_A) {
            count += n_choose_2(GET_A2PATH_ENTRY(g, u, j));
          } else {
            count += n_choose_2(GET_B2PATH_ENTRY(g, u, j));
          }
        }
      }
    }
  } else {
    /* iterate over all nodes that are distance 2 from u */
    for (k = 0; k < g->degree[u]; k++){
      i = g->edgelist[u][k];
      for (l = 0; l < g->degree[i]; l++) {
        j = g->edgelist[i][l];
        if (j != u && !visited[j]) {
          visited[j] = TRUE;
          count += n_choose_2(GET_2PATH_ENTRY(g, u, j));
        }
      }
    }
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
 * Change statistic for Edge
 */
double changeEdge(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)g; (void)i; (void)j; (void)lambda; /* unused parameters */
  assert(!g->is_directed);
  slow_assert(!isEdge(g, i, j));
  return 1;
}

/*
 * Change statistic for 2-stars
 */
double changeTwoStars(graph_t *g, uint_t i, uint_t j, double lambda)
{
  (void)lambda; /* unused parameters */
  assert(!g->is_directed);
  slow_assert(!isEdge(g, i, j));
  return (double)(change_s_stars(g, i, 2) + change_s_stars(g, j, 2));
}

/*
 * Change statistic for alternating k-stars (AS)
 */
double changeAltStars(graph_t *g, uint_t i, uint_t j, double lambda)
{
  assert(lambda > 1);
  assert(!g->is_directed);
  slow_assert(!isEdge(g, i, j));
  return lambda * (2 -
                   POW_LOOKUP(1-1/lambda, g->degree[i]) -
                   POW_LOOKUP(1-1/lambda, g->degree[j]));
}


/*
 * Change statistics for alternating two-path (A2P)
 */
double changeAltTwoPaths(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t v,k;
  double delta = 0;
  assert(lambda > 1);
  assert(!g->is_directed);
  slow_assert(!isEdge(g, i, j));

  if (i == j) {
    return 0;
  }

  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    if (v == i || v == j)
      continue;
    delta += POW_LOOKUP(1-1/lambda, GET_2PATH_ENTRY(g, i, v));
  }
  for (k = 0; k < g->degree[i]; k++) {
    v = g->edgelist[i][k];
    if (v == i || v == j)
      continue;
    delta += POW_LOOKUP(1-1/lambda, GET_2PATH_ENTRY(g, j, v));
  }

  return delta;
}


/*
 * Change statistic for alternating k-triangles (AT)
 */
double changeAltKTriangles(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t v,k,tmp;
  double  delta = 0;
  assert(lambda > 1);
  assert(!g->is_directed);
  slow_assert(!isEdge(g, i, j));
  
  if (i == j) {
    return 0;
  }
  if (g->degree[i] < g->degree[j]) {
    tmp = i;
    i = j;
    j = tmp;
  }

  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    if (v == i || v == j)
      continue;
    if (isEdge(g, i, v))
      delta += POW_LOOKUP(1-1/lambda, GET_2PATH_ENTRY(g, i, v)) +
        POW_LOOKUP(1-1/lambda, GET_2PATH_ENTRY(g, v, j));
  }
  delta += lambda * (1 - POW_LOOKUP(1-1/lambda, GET_2PATH_ENTRY(g, i, j)));
  return delta;
}


/*
 * Change statistic for 4-cycles (4-cycle in PNet, Cycle4A or Cycle4B
 * in MPNet). Note can also be used as for bipartite 4-cycle (C4 in
 * BPNet, X4Cycle in MPNet)
 */
double changeFourCycles(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t v,k,tmp;
  ulong_t delta = 0;
  (void)lambda; /* unused parameters */
  slow_assert(!isEdge(g, i, j));

  if (i == j) {
    return 0;
  }

  /* iterate over neighbours of node with smaller degree */
  if (g->degree[i] < g->degree[j]) {
    tmp = i;
    i = j;
    j = tmp;
  }

  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    if (v == i || v == j) {
      continue;
    }
    if (g->is_bipartite) {
      if (bipartite_node_mode(g, j) == MODE_A) {
	assert(bipartite_node_mode(g, v) == MODE_B);
	delta += GET_B2PATH_ENTRY(g, i, v);
      } else {
	assert(bipartite_node_mode(g, v) == MODE_A);
	delta += GET_A2PATH_ENTRY(g, i, v);
      }
    } else {
      assert(GET_2PATH_ENTRY(g, i, v) == GET_2PATH_ENTRY(g, v, i));
      delta += GET_2PATH_ENTRY(g, i, v);
    }
  }
  
  return (double)delta;
}

/*
 * Change statistic for  3-paths (3-path or D3 in PNet).
 * Note can also be used for bipartite 3-path (L3 in BPNet, X3Path in MPNet)
 */
double changeThreePaths(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t  v,k;
  ulong_t delta = g->degree[i] * g->degree[j];
  (void)lambda; /* unused parameters */
  slow_assert(!isEdge(g, i, j));

  if (i == j) {
    return 0;
  }

  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    delta += g->degree[v] - 1;
  }
  for (k = 0; k < g->degree[i]; k++) {
    v = g->edgelist[i][k];
    delta += g->degree[v] - 1;
  }

  return (double)delta;
}

/*
 * Change statistic for isolated edges.
 * An isolated edge is an edge connecting two nodes each of which has degree 1.
 */
double changeIsolateEdges(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t  v;
  long    delta = 0; /* signed as can be negative */
  (void)lambda; /* unused parameters */
  slow_assert(!isEdge(g, i, j));

  if (i == j) {
    return 0;
  }

  /* if i and j both have degree 0 this created an isolated edge */
  if (g->degree[i] == 0 && g->degree[j] == 0) {
    return 1;
  }

  /* if i has degree 1 and its (therefore only) neighbour v also has degree 1 
   * then adding an edge to i changes the isolated edge i--v to a
   * non-isolated edge. (Similarly for j).
   */
  if (g->degree[i] == 1) {
    v = g->edgelist[i][0];
    if (v != i && v != j && g->degree[v] == 1) {
      delta--;
    }
  }
  if (g->degree[j] == 1) {
    v = g->edgelist[j][0];
    if (v != i && v != j && g->degree[v] == 1) {
      delta--;
    }
    
  }
  return (double)delta;
}

/************************* Actor attribute (binary) **************************/


/*
 * Change statistic for Activity
 */
double changeActivity(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent)
{
  (void)isDelete; /*unused parameters*/
  (void)exponent; /*unused parameters*/
  assert(!g->is_directed);
  slow_assert(!isEdge(g, i, j));
  return ((g->binattr[a][i] == BIN_NA ? 0 : g->binattr[a][i]) +
          (g->binattr[a][j] == BIN_NA ? 0 : g->binattr[a][j]));


}

/******************Attribute interaction (binary) ***********************/

/*
 * Change statistic for interaction effect of two different binary attributes:
 * adding an arc i->j increase the statistic by 1 when binary attribute a
 * on node i is true AND binary attribute b on node j is true.
 *
 * Note that this is particularly useful on bipartite (two-mode) networks,
 * where the two different modes have different attributes. In the bipartite
 * case, attribute a is for mode A nodes, and attribute B is for mode B nodes.
 */
double changeBinaryPairInteraction(graph_t *g, uint_t i, uint_t j,
                                   uint_t a, uint_t b)
{
  uint_t tmp;
  assert(!g->is_directed);
  slow_assert(!isArcOrEdge(g, i, j));
  if (g->is_bipartite) {
    if (bipartite_node_mode(g, i) == MODE_B) {
      assert(bipartite_node_mode(g, j) == MODE_A);
      /* i is mode B node, so swap a and b,
         so that b is used for i and a for node j*/
      tmp = b;
      b = a;
      a = tmp;
    }
  }
  return g->binattr[a][i] != BIN_NA && g->binattr[b][j] != BIN_NA &&
    g->binattr[a][i] && g->binattr[b][j];
}



/*****************************************************************************
 *
 * experimental change statistics functions
 *
 ****************************************************************************/


/*
 * Change statistic for 4-cycles raised to a power. The lambda
 * parameter (> 1.0) (mis)used to specify the value 1/lambda as the
 * epxonent. Note this is not the same meaning of lambda as its
 * original use in the "alternating" parameters.
 *
 * Note can also be used as for bipartite networks.
 *
 */
double changePowerFourCycles(graph_t *g, uint_t i, uint_t j, double lambda)
{
  uint_t  v,k,tmp;
  ulong_t delta = 0;
  ulong_t count = 0, ncount = 0;
  double  alpha = 1/lambda;
  double  change = 0;

  slow_assert(!isEdge(g, i, j));

  if (i == j) {
    return 0;
  }

  /* Number of four-cycles that nodes i or j are already inolved in */
  ulong_t count_i =  num_four_cycles_node(g, i);
  ulong_t count_j =  num_four_cycles_node(g, j);
  count = count_i + count_j;

  /* change statistic for four-cycles */
  delta = changeFourCycles(g, i, j, lambda);
  change = pow(count_i + delta, alpha) - pow(count_i, alpha) +
    pow(count_j + delta, alpha) - pow(count_j, alpha);

  /* neighbours of i */
  for (k = 0; k < g->degree[i]; k++) {
    v = g->edgelist[i][k];
    ncount = num_four_cycles_node(g, v);

    /* TODO compute delta directly instead of counting with/without edge */
    insertEdge(g, i, j);
    uint newcount = num_four_cycles_node(g, v);
    removeEdge(g, i, j);
    change += pow(newcount, alpha) - pow(ncount, alpha);
  }

  /* neighbours of j that are not also neibhours of i and so already counted */
  for (k = 0; k < g->degree[j]; k++) {
    v = g->edgelist[j][k];
    if (!isEdge(g, v, i)) {
      ncount = num_four_cycles_node(g, v);
      
      /* TODO compute delta directly instead of counting with/without edge */
      insertEdge(g, i, j);
      uint newcount = num_four_cycles_node(g, v);
      removeEdge(g, i, j);
      change += pow(newcount, alpha) - pow(ncount, alpha);
    }
  }
  return change;
}
