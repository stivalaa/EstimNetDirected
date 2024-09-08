/*****************************************************************************
 * 
 * File:    bipartiteStats.c
 * Author:  Alex Stivala
 * Created: June 2024
 *
 *
 * statistics functions (summing change statistics is verified against these)
 *
 *
 ****************************************************************************/

#include <assert.h>
#include <math.h>
#include "changeStatisticsUndirected.h"
#include "bipartiteStats.h"

/*****************************************************************************
 *
 * utility functions
 *
 ****************************************************************************/


/*
 * count k-two-paths, as defined by eqn (6.11) in
 *
 *   Wang, P., Sharpe, K., Robins, G. L., & Pattison,
 *   P. E. (2009). Exponential random graph (p∗) models for affiliation
 *   networks. Social Networks, 31(1), 12-25.
 *
 */
static ulonglong_t k_two_paths_A(const graph_t *g, uint_t k)
{
  uint_t l,i;
  ulonglong_t count = 0;

  assert(k > 0);

  for (l = g->num_A_nodes + 1; l < g->num_A_nodes + g->num_B_nodes; l++){
    for (i = g->num_A_nodes; i < l; i++) {
      assert(bipartite_node_mode(g, i) == MODE_B);
      assert(bipartite_node_mode(g, l) == MODE_B);
      count += n_choose_k(GET_B2PATH_ENTRY(g, i, l), k);
    }
  }
  /* Note despite eqn (6.11) in Wang et al. (2009) having this
     statistic multiplied by 1/2 when k = 2 "due to symmetry", this is
     not actually correct for bipartite networks here as we are
     considering only either mode A or mdoe B nodes (not both at
     once). So the symmetry which exists for one-mode networks as per
     Snijders et al. (2006) "New specifications for exponential random
     graph models" [eqns (25a,b) and (26a), pp. 123-124] when the
     summation over all i < j means two pairs of nodes are considered
     in a four-cycle, is not true here as the summation over i < l
     considers only one pair of nodes (mode B nodes in k_two_paths_A()
     or mode A nodes in k_two_paths_B) are considered.
  */
  return count;
}

static ulonglong_t k_two_paths_B(const graph_t *g, uint_t k)
{
  uint_t l,i;
  ulonglong_t count = 0;

  assert(k > 0);

  for (l = 1; l < g->num_A_nodes; l++){
    for (i = 0; i < l; i++) {
      assert(bipartite_node_mode(g, i) == MODE_A);
      assert(bipartite_node_mode(g, l) == MODE_A);
      count += n_choose_k(GET_A2PATH_ENTRY(g, i, l), k);
    }
  }
  /* Note despite eqn (6.11) in Wang et al. (2009) having this
     statistic multiplied by 1/2 when k = 2 "due to symmetry", this is
     not actually correct for bipartite networks here as we are
     considering only either mode A or mdoe B nodes (not both at
     once). So the symmetry which exists for one-mode networks as per
     Snijders et al. (2006) "New specifications for exponential random
     graph models" [eqns (25a,b) and (26a), pp. 123-124] when the
     summation over all i < j means two pairs of nodes are considered
     in a four-cycle, is not true here as the summation over i < l
     considers only one pair of nodes (mode B nodes in k_two_paths_A()
     or mode A nodes in k_two_paths_B) are considered.
  */
  return count;
}


/*
 * Count number of four-cycles that a particular node u is involved in.
 *
 * This version for bipartite networks.
 */
static uint_t num_four_cycles_node_SLOW(const graph_t *g, uint_t u)
{
  uint_t v;
  uint_t count = 0;

  if (g->is_bipartite) {
    if (bipartite_node_mode(g, u) == MODE_A) {
      for (v = 0; v < g->num_A_nodes; v++) {
        if (v != u) {
          assert(bipartite_node_mode(g, v) == MODE_A);
          count += n_choose_2(GET_A2PATH_ENTRY(g, u, v));
        }
      }
    } else {
      for (v = g->num_A_nodes; v < g->num_A_nodes + g->num_B_nodes; v++) {
        if (v != u) {
          assert(bipartite_node_mode(g, v) == MODE_B);
          count += n_choose_2(GET_B2PATH_ENTRY(g, u, v));
        }
      }
    }
  }
  return count;
}



/*****************************************************************************
 *
 * statistics functions (summing change statistics is verified against these)
 *
 ****************************************************************************/


/*
 * Statistic for FourCycles, number of four-cycles in a bipartite graph.
 *
 * This version counting over pairs of mode A nodes, but result must be
 * equal to that counting over pairs of mode B nodes instead.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      number of four-cycles in bipartite graph g
 */
double FourCyclesA(const graph_t *g)
{
  uint_t i,l;
  double value = 0;

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = 1; i < g->num_A_nodes; i++) {
    for (l = 0; l < i; l++) {
      assert(bipartite_node_mode(g, i) == MODE_A);
      assert(bipartite_node_mode(g, l) == MODE_A);
      value += n_choose_2(GET_A2PATH_ENTRY(g, i, l));
    }
  }
  return value;
}

/*
 * Statistic for FourCycles, number of four-cycles in a bipartite graph.
 *
 * This version counting over pairs of mode B nodes, but result must be
 * equal to that counting over pairs of mode A nodes instead.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      number of four-cycles in bipartite graph g
 */
double FourCyclesB(const graph_t *g)
{
  uint_t i,l;
  double value = 0;

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = g->num_A_nodes + 1; i < g->num_A_nodes + g->num_B_nodes; i++) {
    for (l = g->num_A_nodes; l < i; l++) {
      assert(bipartite_node_mode(g, i) == MODE_B);
      assert(bipartite_node_mode(g, l) == MODE_B);
      value += n_choose_2(GET_B2PATH_ENTRY(g, i, l));
    }
  }
  return value;
}


/*
 * Statistic for BipartiteAltKCyclesA, alternating k-cycles for type A
 * nodes (K-Ca in BPNet, XACA in MPNet) defined by eqn (6.12) in:
 *
 *   Wang, P., Sharpe, K., Robins, G. L., & Pattison,
 *   P. E. (2009). Exponential random graph (p∗) models for affiliation
 *   networks. Social Networks, 31(1), 12-25.
 *
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *     lambda - decay value > 1.0
 *
 * Return value:
 *      statistic for g with decay value lambda
 */
double BipartiteAltKCyclesA(const graph_t *g, double lambda)
{
  uint_t i,l;
  double value = 0;

  assert (lambda > 1.0);

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (l = g->num_A_nodes + 1; l < g->num_A_nodes + g->num_B_nodes; l++) {
    for (i = g->num_A_nodes; i < l; i++) {
      assert(bipartite_node_mode(g, i) == MODE_B);
      assert(bipartite_node_mode(g, l) == MODE_B);
      value += 1 - POW_LOOKUP(1-1/lambda, GET_B2PATH_ENTRY(g, i, l));
    }
  }
  return lambda * value;
}


/*
 * Statistic for BipartiteAltKCyclesB, alternating k-cycles for type B
 * nodes (K-Cp in BPNet, XACB in MPNet) defined by eqn (6.12) in:
 *
 *   Wang, P., Sharpe, K., Robins, G. L., & Pattison,
 *   P. E. (2009). Exponential random graph (p∗) models for affiliation
 *   networks. Social Networks, 31(1), 12-25.
 *
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *     lambda - decay value > 1.0
 *
 * Return value:
 *      statistic for g with decay value lambda
 */
double BipartiteAltKCyclesB(const graph_t *g, double lambda)
{
  uint_t i,l;
  double value = 0;

  assert (lambda > 1.0);

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (l = 1; l < g->num_A_nodes; l++) {
    for (i = 0; i < l; i++) {
      assert(bipartite_node_mode(g, i) == MODE_A);
      assert(bipartite_node_mode(g, l) == MODE_A);
      value += 1 - POW_LOOKUP(1-1/lambda, GET_A2PATH_ENTRY(g, i, l));
    }
  }
  return lambda * value;
}




/*
 * Statistic for BipartiteAltKCyclesA, alternating k-cycles for type A
 * nodes (K-Ca in BPNet, XACA in MPNet), alternative (inefficient)
 * implementation, defined by eqn (6.12) in:
 *
 *   Wang, P., Sharpe, K., Robins, G. L., & Pattison,
 *   P. E. (2009). Exponential random graph (p∗) models for affiliation
 *   networks. Social Networks, 31(1), 12-25.
 *
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *     lambda - decay value > 1.0
 *
 * Return value:
 *      statistic for g with decay value lambda
 */
double BipartiteAltKCyclesA_SLOW(const graph_t *g, double lambda)
{
  uint_t i;
  double value;

  assert (lambda > 1.0);

  assert(g->is_bipartite);
  assert(!g->is_directed);

  /* Note despite eqn (6.12) in Wang et al. (2009) multipliying the
     second term [k_two_paths_A(g, 2)/lambda] by 2, this is not
     actually correct as the division by two in eqn (6.11) is not
     correct (see comment in k_two_paths_A()), so there is no factor
     of 2 here.
  */
  value = k_two_paths_A(g, 1) - k_two_paths_A(g, 2)/lambda;

  for (i = 3; i < g->num_A_nodes + g->num_B_nodes - 1; i++) {
    value += pow(-1/lambda, i-1) * k_two_paths_A(g, i);
  }
  return value;
}

/*
 * Statistic for BipartiteAltKCyclesB, alternating k-cycles for type B
 * nodes (K-Cp in BPNet, XACB in MPNet), alternative (inefficient)
 * implementation, defined by eqn (6.12) in:
 *
 *   Wang, P., Sharpe, K., Robins, G. L., & Pattison,
 *   P. E. (2009). Exponential random graph (p∗) models for affiliation
 *   networks. Social Networks, 31(1), 12-25.
 *
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *     lambda - decay value > 1.0
 *
 * Return value:
 *      statistic for g with decay value lambda
 */
double BipartiteAltKCyclesB_SLOW(const graph_t *g, double lambda)
{
  uint_t i;
  double value;

  assert (lambda > 1.0);

  assert(g->is_bipartite);
  assert(!g->is_directed);

  /* Note despite eqn (6.12) in Wang et al. (2009) multipliying the
     second term [k_two_paths_A(g, 2)/lambda] by 2, this is not
     actually correct as the division by two in eqn (6.11) is not
     correct (see comment in k_two_paths_B()), so there is no factor
     of 2 here.
  */
    value = k_two_paths_B(g, 1) - k_two_paths_B(g, 2)/lambda;

  for (i = 3; i < g->num_A_nodes + g->num_B_nodes - 1; i++) {
    value += pow(-1/lambda, i-1) * k_two_paths_B(g, i);
  }
  return value;
}



/*
 * Statistic for alternating k-4-cycles for type A nodes (new change
 * statistic suggested in email (basically paper outline, with
 * spreadsheet attachments for literature search, examples, etc.)
 * "idea for a (slightly) new bipartite change statistic" sent 23 Nov
 * 2022):
 *
 *   The proposed new statistic is a very simple modification of the
 *   "Alternating k-two-paths" bipartite statistic (K-Ca and K-Cp)
 *   statistics described by Wang et al. (2009, p.19). I propose to
 *   simply remove the first term of Wang et al. (2009) equation 6.12,
 *   and reverse the signs, so that it no longer counts open two-paths,
 *   but the first, positive, term actually counts four-cycles.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *     lambda - decay value > 1.0
 *
 * Return value:
 *      statistic for g with decay value lambda
 */
double BipartiteAltK4CyclesA_SLOW(const graph_t *g, double lambda)
{
  uint_t i;
  double value;

  assert (lambda > 1.0);

  assert(g->is_bipartite);
  assert(!g->is_directed);

  value = k_two_paths_A(g, 2)/lambda;

  for (i = 3; i < g->num_A_nodes + g->num_B_nodes - 1; i++) {
    value += -1 * pow(-1/lambda, i-1) * k_two_paths_A(g, i);
  }
  return value;
}



/*
 * Statistic for number of 4-cycles at each node raised to a
 * power. The lambda parameter (> 1.0) (mis)used to specify the value
 * 1/lambda as the epxonent. Note this is not the same meaning of
 * lambda as its original use in the "alternating" parameters.
 *
 * This version counting over pairs of mode A nodes only.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      
 */
double PowerFourCyclesA(const graph_t *g, double lambda)
{
  uint_t  i;
  double  alpha = 1/lambda;
  double  value = 0;
  uint_t  fourcycle_count = 0;
  ulonglong_t fourcycle_count_sum = 0;
  ulonglong_t num_fourcycles;

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = 0; i < g->num_A_nodes; i++) {
    uint_t fourcycle_count_SLOW = num_four_cycles_node_SLOW(g, i);
    fourcycle_count = num_four_cycles_node(g, i);
    assert(fourcycle_count == fourcycle_count_SLOW);
    fourcycle_count_sum += fourcycle_count;
    value += pow(num_four_cycles_node(g, i), alpha);
  }
  /* Each four-cycle contains 2 nodes in mode A so is counted twice */
  assert(fourcycle_count_sum % 2 == 0);
  num_fourcycles = FourCyclesA(g);
  assert(fourcycle_count_sum / 2 == num_fourcycles);
  return value;
}


/*
 * Statistic for number of 4-cycles at each node raised to a
 * power. The lambda parameter (> 1.0) (mis)used to specify the value
 * 1/lambda as the epxonent. Note this is not the same meaning of
 * lambda as its original use in the "alternating" parameters.
 *
 * This version counting over pairs of mode B nodes only.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      
 */
double PowerFourCyclesB(const graph_t *g, double lambda)
{
  uint_t  i;
  double  alpha = 1/lambda;
  double  value = 0;
  uint_t  fourcycle_count = 0;
  ulonglong_t fourcycle_count_sum = 0;
  ulonglong_t num_fourcycles;

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = g->num_A_nodes; i < g->num_A_nodes + g->num_B_nodes; i++) {
    uint_t fourcycle_count_SLOW = num_four_cycles_node_SLOW(g, i);
    fourcycle_count = num_four_cycles_node(g, i);
    assert(fourcycle_count == fourcycle_count_SLOW);
    fourcycle_count_sum += fourcycle_count;
    value += pow(num_four_cycles_node(g, i), alpha);
  }
  /* Each four-cycle contains 2 nodes in mode B so is counted twice */
  assert(fourcycle_count_sum % 2 == 0);
  num_fourcycles = FourCyclesB(g);
  assert(fourcycle_count_sum / 2 == num_fourcycles);
  return value;
}



/*
 * Statistic for count of 4-cycles at each node raised to a power. The
 * lambda parameter (> 1.0) (mis)used to specify the value 1/lambda as
 * the epxonent. Note this is not the same meaning of lambda as its
 * original use in the "alternating" parameters.
 *
 * Parameters:
 *     g      - undirected bipartite graph
 *
 * Return value:
 *      number of four-cycles in bipartite graph g
 */
double PowerFourCycles(const graph_t *g, double lambda)
{
  assert(g->is_bipartite);
  assert(!g->is_directed);

  double power_4cycles_A = PowerFourCyclesA(g, lambda);
  double power_4cycles_B = PowerFourCyclesB(g, lambda);
  return power_4cycles_A + power_4cycles_B;
}

