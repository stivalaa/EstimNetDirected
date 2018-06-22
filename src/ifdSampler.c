/*****************************************************************************
 * 
 * File:    ifdSampler.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: June 2018
 *
 * Improved fixed density (IFD) ERGM distribution sampler
 *
 * Byshkin, M., Stivala, A., Mira, A., Krause, R., Robins, G., & Lomi,
 * A. (2016). Auxiliary parameter MCMC for exponential random graph
 * models. Journal of Statistical Physics, 165(4), 740-754.
 * the arc.
 *
 *
 ****************************************************************************/

#include <assert.h>
#include <math.h>
#include "utils.h"
#include "changeStatisticsDirected.h"
#include "ifdSampler.h"

/*
 * Return the value to subtract from the IFD auxiliary parameter in
 * order to get the Arc parameter value when using the IFD sampler,
 * i.e. recover the parameter Theta_L when we have the value of the
 * auxiliary parameter V from equation (20) in the paper (but for
 * directed graphs here).
 *
 * Parameters:
 *    g - digraph object
 * 
 * Return value:
 *    edge correction value
 */
double arcCorrection(const digraph_t *g) {
  double N         = g->num_nodes;
  double num_dyads = N*(N-1);/*directed so not div by 2*/
  double num_arcs  = g->num_arcs;
  return log((num_dyads - num_arcs) / (num_arcs + 1));
}


/*
 * Improved Fixed Density (IFD) ERGM MCMC sampler, described in:
 * 
 * Byshkin, M., Stivala, A., Mira, A., Krause, R., Robins, G., & Lomi,
 * A. (2016). Auxiliary parameter MCMC for exponential random graph
 * models. Journal of Statistical Physics, 165(4), 740-754.
 *
 *
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
 *   ifd_K       - constant for multipliying step size of auxiliary parameter
 *   dzArc       - (Out) Arc statistic differnce from observed: just Ndel-Nadd
 *   ifd_aux_param  - (In/Out) IFD auxiliary parameter. Pass zero initially, then
 *                    reuse each call to update.
 *
 * Return value:
 *   Acceptance rate.
 *
 * The addChangeStats and delChangeStats arrays are of length n corresponding
 * to the theta parameter array and change_stats_funcs change statistics
 * function pointer array. On exit they are set to the sum values of the
 * change statistics for add and delete moves respectively.
 */
double ifdSampler(digraph_t *g,  uint_t n, uint_t n_attr, uint_t n_dyadic,
                  change_stats_func_t *change_stats_funcs[],
                  attr_change_stats_func_t *attr_change_stats_funcs[],
                  dyadic_change_stats_func_t *dyadic_change_stats_funcs[],
                  uint_t attr_indices[], double theta[],
                  double addChangeStats[], double delChangeStats[],
                  uint_t sampler_m,
                  bool performMove,
                  double ifd_K, double *dzArc, double *ifd_aux_param)
{
  static bool   isDelete = FALSE; /* delete or add move. FIXME don't use static, make param */

  double *changestats = (double *)safe_malloc(n*sizeof(double));
  double  total;        /* sum of theta*changestats */
  uint_t  accepted = 0; /* number of accepted moves */
  /* Ndel and Nadd are int not uint_t as we do signed math with them */
  int     Ndel = 0;     /* number of add moves */
  int     Nadd = 0;     /* number of delete moves */
  double  ifd_aux_param_step;
  double  acceptance_rate;
  uint_t  i,j,k,l,param_i;
  uint_t  arcidx = 0;

  for (i = 0; i < n; i++)
    addChangeStats[i] = delChangeStats[i] = 0;

  for (k = 0; k < sampler_m; k++) {
    if (isDelete) {
      /* Delete move. Find an exising arc uniformly at random to delete. */
      arcidx = int_urand(g->num_arcs);
      i = g->allarcs[arcidx].i;
      j = g->allarcs[arcidx].j;
      assert(isArc(g, i, j));
    } else {
      /* Add move. Find two nodes i, j without arc i->j uniformly at
         random. Because graph is sparse, it is not too inefficient
         to just pick random nodes until such a pair is found */
      do {
        i = int_urand(g->num_nodes);
        do {
          j = int_urand(g->num_nodes);
        } while (i == j);
      } while (isArc(g, i, j));
    }
    
    /* The change statistics are all computed on the basis of adding arc i->j
       so for a delete move, we (perhaps temporarily) remove it to compute the
       change statistics, and negate them */
    SAMPLER_DEBUG_PRINT(("%s %d -> %d\n",isDelete ? "del" : "add", i, j));
    if (isDelete) {
      removeArc_allarcs(g, i, j, arcidx);
      Ndel++;
    } else {
      Nadd++;
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
    
    /* add the IFD auxiliary parameter value */
    total += (isDelete ? -1 : 1) * *ifd_aux_param;

    /* now exp(total) is the acceptance probability */
    if (urand() < exp(total)) {
      accepted++;
      if (performMove) {
        /* actually do the move. If deleting, already done it. For add, add
           the arc now */
        if (!isDelete) {
          assert(!isArc(g, i, j));
          insertArc_allarcs(g, i, j);
        }
      } else {
        /* not actually doing the moves, so reverse change for delete move
           to restore g to original state */
        if (isDelete) {
          assert(!isArc(g, i, j));
          insertArc_allarcs(g, i, j);
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
      isDelete = !isDelete;
    } else {
      /* move not acceptd, so reverse change for delete */
      if (isDelete) {
        assert(!isArc(g, i, j));
        insertArc_allarcs(g, i, j);
      }
    }
  }
  
  /* update the IFD auxiliary parameter */
  ifd_aux_param_step = ifd_K * (Ndel - Nadd) * (Ndel - Nadd) / ((Ndel+Nadd)*(Ndel+Nadd));
  SAMPLER_DEBUG_PRINT(("accepted = %u, Nadd = %u, Ndel = %u, ifd_aux_param_step = %g\n",
                       accepted, Nadd, Ndel, ifd_aux_param_step));
  if (Ndel - Nadd > 0) {
    *ifd_aux_param -= ifd_aux_param_step;
  } else if (Ndel - Nadd < 0) {
    *ifd_aux_param += ifd_aux_param_step;
  }
  SAMPLER_DEBUG_PRINT(("ifd_aux_param = %g\n", *ifd_aux_param));
  if (fabs(Ndel - Nadd) / (Ndel + Nadd) > 0.8) {
    fprintf(stderr,
            "WARNING: IFD sampler Ndel = %d Nadd = %d ifd_aux_param = %g increase ifd_K = %f\n",
            Ndel, Nadd, *ifd_aux_param, ifd_K);
  }

  *dzArc = (double)Ndel - (double)Nadd;
  acceptance_rate = (double)accepted / sampler_m;
  free(changestats);
  return acceptance_rate;
}
