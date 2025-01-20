/*****************************************************************************
 * 
 * File:    attrbipartiteStats.c
 * Author:  Alex Stivala
 * Created: January 2025
 *
 *
 * statistics functions (summing change statistics is verified against these)
 *
 *
 ****************************************************************************/

#include <assert.h>
#include <math.h>
#include "changeStatisticsBipartiteUndirected.h"
#include "attrbipartiteStats.h"


/*****************************************************************************
 *
 * utility functions
 *
 ****************************************************************************/

/*
 * Return the number of neighbours of node i with binary attribute a.
 */
static uint_t count_neighbours_with_binattr_a(const graph_t *g, uint_t i, uint_t a)
{
  uint_t num_neighbours_with_a = 0;
  uint_t k, v;

  assert(!g->is_directed);

  for (k = 0; k < g->degree[i]; k++) {
    v = g->edgelist[i][k];
    if (g->binattr[a][v] != BIN_NA && g->binattr[a][v]) {
      num_neighbours_with_a++;
    }
  }
  return num_neighbours_with_a;
}


/*****************************************************************************
 *
 * statistics functions (summing change statistics is verified against these)
 *
 ****************************************************************************/

/*
 * Statistic for bipartite exactly one neighbour with binary attribute a
 * for type A nodes.
 *
 * The statistic counts the number of type A nodes that have exactly one
 * neighbour (therefore of type B) with the binary attribute a.
 *
 * Note that binary attribute a here is a binary attribute for type B nodes.
 */
double BipartiteExactlyOneNeighbourA(const graph_t *g, uint_t a)

{
  uint_t i;
  uint_t value = 0;

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = 0; i < g->num_A_nodes; i++) {
    assert(bipartite_node_mode(g, i) == MODE_A);
    value += count_neighbours_with_binattr_a(g, i, a) == 1 ? 1 : 0;
  }
  return (double)value;
}


/*
 * Statistic for bipartite exactly one neighbour with binary attribute a
 * for type B nodes.
 *
 * The statistic counts the number of type B nodes that have exactly one
 * neighbour (therefore of type A) with the binary attribute a.
 *
 * Note that binary attribute a here is a binary attribute for type A nodes.
 */
double BipartiteExactlyOneNeighbourB(const graph_t *g, uint_t a)

{
  uint_t i;
  uint_t value = 0;

  assert(g->is_bipartite);
  assert(!g->is_directed);

  for (i = g->num_A_nodes; i < g->num_nodes; i++) {
    assert(bipartite_node_mode(g, i) == MODE_B);
    value += count_neighbours_with_binattr_a(g, i, a) == 1 ? 1 : 0;
  }
  return (double)value;
}


/*
 * Statistic for Bipartite 2-path beween two type A nodes
 * each of which has  exactly one neighbour with binary attribute a.
 *
 * The statistic counts the number of two-paths between pairs of type
 * A nodes that both have exactly one neighbour (therefore of type B) with
 * the binary attribute a.
 *
 * Note that binary attribute a here is a binary attribute for type B nodes.
 */
double BipartiteTwoPathExactlyOneNeighbourA(const graph_t *g, uint_t a)
{
  uint_t i,j,k,l,v;
  uint_t count = 0;

  assert(g->is_bipartite);
  assert(!g->is_directed);
  for (i = 0; i < g->num_A_nodes; i++) {
    assert(bipartite_node_mode(g, i) == MODE_A);
    if (count_neighbours_with_binattr_a(g, i, a) != 1)
      continue;
    for (j = i+1; j < g->num_A_nodes; j++) {
      assert(bipartite_node_mode(g, j) == MODE_A);
      if (count_neighbours_with_binattr_a(g, j, a) != 1)
        continue;
      for (k = 0; k < g->degree[i]; k++)  {
        v = g->edgelist[i][k];   /* i -- v */
        assert(bipartite_node_mode(g, v) == MODE_B);
        for (l = 0; l < g->degree[j]; l++) {
          if (g->edgelist[j][l] == v) {   /* v -- j */
            count++;
          }
        }
      }
    }
  }
  return (double)count;
}
