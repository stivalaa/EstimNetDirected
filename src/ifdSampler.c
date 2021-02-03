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
 *
 * It also optionally does conditional estimation for snowball sampled
 * network. In this case in the MCMC algorithm the ties between nodes
 * in the outermost wave are fixed, as are ties between nodes in the
 * outermost wave and the preceding (second-last) wave. In addition, a
 * tie cannot be added if it would "skip over" a wave (i.e. the
 * absolute difference in wave number between the nodes to add a tie must
 * be at most 1), and a tie cannot be deleted if it is the last remaining
 * tie connected a node to the preceding wave.
 *
 * Note in the case of directed networks the snowball sampling procedure
 * has been assumed to ignore the direction of arcs, so when we consider
 * the above rules here we ignore the direction of the arcs also.
 *
 * References for conditional estimation of snowball sampled network are:
 * 
 * Pattison, P. E., Robins, G. L., Snijders, T. A., & Wang,
 * P. (2013). Conditional estimation of exponential random graph
 * models from snowball sampling designs. Journal of Mathematical
 * Psychology, 57(6), 284-296.
 * 
 * Stivala, A. D., Koskinen, J. H., Rolls, D. A., Wang, P., & Robins,
 * G. L. (2016). Snowball sampling for estimating exponential random
 * graph models for large networks. Social Networks, 47, 167-188.
 *
 * And for the directed networks case specifically:
 *
 * Stivala, A., Rolls, D., & Robins, G. (2015). The ins and outs of
 * snowball sampling: ERGM estimation for very large directed
 * networks, presented at INSNA Sunbelt XXXV Conference, Brighton UK,
 * June 23-28, 2015. [Slides available from
 * https://sites.google.com/site/alexdstivala/home/conferences]
 *
 * Stivala, A., Rolls, D., & Robins, G. (2018). Estimating exponential
 * random graph models for large directed networks with snowball
 * sampling. Unpublished manuscript.
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
 *   n_attr_interaction - number of attribute interaction change stats funcs
 *   change_stats_funcs - array of pointers to change statistics functions
 *                        length is n-n_attr-n_dyadic-n_attr_interaction
 *   attr_change_stats_funcs - array of pointers to change statistics functions
 *                             length is n_attr
 *   dyadic_change_stats_funcs - array of pointers to dyadic change stats funcs
 *                             length is n_dyadic
 *   attr_interaction_change_stats_funcs - array of pointers to attribute
 *                           interaction (pair) change statistics functions.
 *                           length is n_attr_interaction.
 *   attr_indices   - array of n_attr attribute indices (index into g->binattr
 *                    or g->catattr) corresponding to attr_change_stats_funcs
 *                    E.g. for Sender effect on the first binary attribute,
 *                    attr_indices[x] = 0 and attr_change_stats_funcs[x] =
 *                    changeSender
 *   attr_interaction_pair_indices - array of n_attr_interaction pairs
 *                          of attribute inidices similar to above but
 *                          for attr_interaction_change_setats_funcs which
 *                          requires pairs of indices.
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
 *   useConditionalEstimation - if True do conditional estimation of snowball
 *                              network sample.
 *   forbidReciprocity - if True do not allow reciprocated arcs.
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
                  uint_t n_attr_interaction,
                  change_stats_func_t *change_stats_funcs[],
                  attr_change_stats_func_t *attr_change_stats_funcs[],
                  dyadic_change_stats_func_t *dyadic_change_stats_funcs[],
                  attr_interaction_change_stats_func_t
                                   *attr_interaction_change_stats_funcs[],
                  uint_t attr_indices[],
                  uint_pair_t attr_interaction_pair_indices[],
                  double theta[],
                  double addChangeStats[], double delChangeStats[],
                  uint_t sampler_m,
                  bool performMove,
                  double ifd_K, double *dzArc, double *ifd_aux_param,
                  bool useConditionalEstimation,
                  bool forbidReciprocity)
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
  uint_t  i,j,k,l;
  uint_t  arcidx = 0;

  
  for (i = 0; i < n; i++) {
    addChangeStats[i] = delChangeStats[i] = 0;
  }

  for (k = 0; k < sampler_m; k++) {

    if (useConditionalEstimation) {
      assert(!forbidReciprocity); /* TODO not implemented for snowball */
      if (isDelete) {
        /* Delete move for conditional estimation. Find an existing
           arc between nodes in inner waves (i.e. fixing ties in
           outermost wave and between outermost and second-outermost
           waves) uniformly at random to delete.  Extra constraint for
           conditional estimation that a tie cannot be deleted if it
           is last remaining tie connecting node to preceding wave.
           Note ignoring arc direction here as assumed snowball sample
           ignored arc directions.
         */
        do {
          arcidx = int_urand(g->num_inner_arcs);
          i = g->allinnerarcs[arcidx].i;
          j = g->allinnerarcs[arcidx].j;
          SAMPLER_DEBUG_PRINT(("conditional del arcidx %u (%u -> %u) zones %u %u\n", arcidx, i, j, g->zone[i], g->zone[j]));
          assert(g->zone[i] < g->max_zone && g->zone[j] < g->max_zone);
          /* any tie must be within same zone or between adjacent zones */
          assert(labs((long)g->zone[i] - (long)g->zone[j]) <= 1);
        } while ((g->zone[i] > g->zone[j] && g->prev_wave_degree[i] == 1) ||
                 (g->zone[j] > g->zone[i] && g->prev_wave_degree[j] == 1));
      } else {
        /* Add move for conditional estimation. Find two nodes i, j in
           inner waves without arc i->j uniformly at random. Because
           graph is sparse, it is not too inefficient to just pick
           random nodes until such a pair is found. For conditional
           estimation we also have the extra constraint that the nodes
           must be in the same wave or adjacent waves for the tie to
           be added. */
        do {
          i = g->inner_nodes[int_urand(g->num_inner_nodes)];          
          do {
            j = g->inner_nodes[int_urand(g->num_inner_nodes)];        
          } while (i == j);
          assert(g->zone[i] < g->max_zone && g->zone[j] < g->max_zone);
        } while (isArc(g, i, j) ||
                 (labs((long)g->zone[i] - (long)g->zone[j]) > 1));
      }
    } else {
      /* not using conditional estimation */
      if (isDelete) {
        /* Delete move. Find an existing arc uniformly at random to delete. */
        arcidx = int_urand(g->num_arcs);
        i = g->allarcs[arcidx].i;
        j = g->allarcs[arcidx].j;
        /*removed as slows significantly: assert(isArc(g, i, j));*/
        /* no need to condsider forbidReciprocity on delete move */
      } else {
        /* Add move. Find two nodes i, j without arc i->j uniformly at
           random. Because graph is sparse, it is not too inefficient
           to just pick random nodes until such a pair is found */
        do {
          do {
            i = int_urand(g->num_nodes);
            do {
              j = int_urand(g->num_nodes);
            } while (i == j);
          } while (isArc(g, i, j));
        } while (forbidReciprocity && isArc(g, j, i));
      }
    }
    
    /* The change statistics are all computed on the basis of adding arc i->j
       so for a delete move, we (perhaps temporarily) remove it to compute the
       change statistics, and negate them */
    SAMPLER_DEBUG_PRINT(("%s %d -> %d\n",isDelete ? "del" : "add", i, j));
    if (isDelete) {
      if (useConditionalEstimation){
        removeArc_allinnerarcs(g, i, j, arcidx);
      } else {
        removeArc_allarcs(g, i, j, arcidx);        
      }
      Ndel++;
    } else {
      Nadd++;
    }

    total = calcChangeStats(g, i, j, n, n_attr, n_dyadic, n_attr_interaction,
                            change_stats_funcs,
                            attr_change_stats_funcs,
                            dyadic_change_stats_funcs,
                            attr_interaction_change_stats_funcs,
                            attr_indices,
                            attr_interaction_pair_indices,
                            theta, isDelete, changestats);

    
    /* add the IFD auxiliary parameter value */
    total += (isDelete ? -1 : 1) * *ifd_aux_param;

    /* now exp(total) is the acceptance probability */
    if (urand() < exp(total)) {
      accepted++;
      if (performMove) {
        /* actually do the move. If deleting, already done it. For add, add
           the arc now */
        if (!isDelete) {
          if (useConditionalEstimation) {
            insertArc_allinnerarcs(g, i, j);
          } else {
            insertArc_allarcs(g, i, j);
          }
        }
      } else {
        /* not actually doing the moves, so reverse change for delete move
           to restore g to original state */
        if (isDelete) {
          if (useConditionalEstimation) {
            insertArc_allinnerarcs(g, i, j);
          } else {
            insertArc_allarcs(g, i, j);
          }
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
        if (useConditionalEstimation) {
          insertArc_allinnerarcs(g, i, j);
        } else {
          insertArc_allarcs(g, i, j);          
        }
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
  if ((double)abs(Ndel - Nadd) / (double)(Ndel + Nadd) > 0.8) {
    fprintf(stderr,
            "WARNING: IFD sampler Ndel = %d Nadd = %d ifd_aux_param = %g increase ifd_K = %f\n",
            Ndel, Nadd, *ifd_aux_param, ifd_K);
  }

  *dzArc = (double)Ndel - (double)Nadd;
  acceptance_rate = (double)accepted / sampler_m;
  free(changestats);
  return acceptance_rate;
}
