/*****************************************************************************
 * 
 * File:    basicSampler.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * "Basic" ERGM distribution sampler. It picks a random dyad and toggles
 * the arc.
 *
 *
 ****************************************************************************/

#include <assert.h>
#include <math.h>
#include "utils.h"
#include "changeStatisticsDirected.h"
#include "basicSampler.h"

/*
 * Basic ERGM MCMC sampler. Uniformly at random a dyad i, j is chosen
 * and the arc i->j is toggled, ie added if it does not exist, removed
 * if it does.
 *
 * Parameters:
 *   g      - digraph object. Modifed if performMove is true.
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
 *   addChangeStats - (Out) vector of n change stats for add moves
 *                    Allocated by caller.
 *   delChangeStats - (Out) vector of n change stats for delete moves
 *                    Allocated by caller
 *   sampler_m   - Number of proposals (sampling iterations)
 *   performMove - if true, moves are actually performed (digraph updated).
 *                 Otherwise digraph is not actually changed.
 *   useConditionalEstimation - if True do conditional estimation of snowball
 *                              network sample.
 *
 * Return value:
 *   Acceptance rate.
 *
 * The addChangeStats and delChangeStats arrays are of length n corresponding
 * to the theta parameter array and change_stats_funcs change statistics
 * function pointer array. On exit they are set to the sum values of the
 * change statistics for add and delete moves respectively.
 *
 * Note that this sampler does not update the digraph allarcs flat arc list
 * as it does not need to use it at all, so that it remains as it was 
 * and therefore becomes inconsistent with the actual graph when it is modified
 * in this function, so should therefore not be used afterwards.
 */
double basicSampler(digraph_t *g,  uint_t n, uint_t n_attr, uint_t n_dyadic,
                    change_stats_func_t *change_stats_funcs[],
                    attr_change_stats_func_t *attr_change_stats_funcs[],
                    dyadic_change_stats_func_t *dyadic_change_stats_funcs[],
                    uint_t attr_indices[], double theta[],
                    double addChangeStats[], double delChangeStats[],
                    uint_t sampler_m,
                    bool performMove,
                    bool useConditionalEstimation)
{
  uint_t accepted = 0;    /* number of accepted moves */
  double acceptance_rate;
  uint_t i,j,k,l,param_i;
  bool   isDelete;
  double *changestats = (double *)safe_malloc(n*sizeof(double));
  double total;  /* sum of theta*changestats */

  for (i = 0; i < n; i++)
    addChangeStats[i] = delChangeStats[i] = 0;

  for (k = 0; k < sampler_m; k++) {
    /* select two nodes i and j uniformly at random and toggle arc between them*/
    i = int_urand(g->num_nodes);
    do {
      j = int_urand(g->num_nodes);
    } while (i == j);


    /* The change statistics are all computed on the basis of adding arc i->j
       so if if the arc exists, we temporarily remove it to compute the
       change statistics, and negate them */
    isDelete = isArc(g, i ,j);
    SAMPLER_DEBUG_PRINT(("%s %d -> %d\n",isDelete ? "del" : "add", i, j));
    if (isDelete) {
      removeArc(g, i, j);
      assert(!isArc(g, i, j));
    }

    total = 0;
    param_i = 0;
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
    

    /* now exp(total) is the acceptance probability */
    if (urand() < exp(total)) {
      accepted++;
      if (performMove) {
        /* actually do the move. If deleting, already done it. For add, add
           the arc now */
        if (!isDelete)
          insertArc(g, i, j);
      } else {
        /* not actually doing the moves, so reverse change for delete move
           to restore g to original state */
        if (isDelete) {
          assert(!isArc(g, i, j));
          insertArc(g, i, j);
        }
      }
      /* accumulate the change statistics for add and del moves separately */
      if (isDelete) {
        for (l = 0; l < n; l++)
          delChangeStats[l] += changestats[l];
      } else {
        for (l = 0; l < n; l++)
          addChangeStats[l] += changestats[l];
      }
    } else {
      /* move not acceptd, so reverse change for delete */
      if (isDelete) {
        assert(!isArc(g, i, j));
        insertArc(g, i, j);
      }
    }
  }
  
  acceptance_rate = (double)accepted / sampler_m;
  free(changestats);
  return acceptance_rate;
}
