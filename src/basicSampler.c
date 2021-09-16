/*****************************************************************************
 * 
 * File:    basicSampler.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * "Basic" ERGM distribution sampler. It picks a random dyad and toggles
 * the arc.
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
 * Also can optionally do citation ERGM (cERGM) estimation, which is 
 * conditional on the term (time period) of the node. All ties
 * except those from a node in the last time period are fixed.
 * 
 * Reference for citation ERGM (cERGM) estimation is:
 *
 *   Schmid, C. S., Chen, T. H. Y., & Desmarais, B. A. (2021). 
 *   Generative Dynamics of Supreme Court Citations:
 *   Analysis with a New Statistical Network Model. arXiv preprint
 *   arXiv:2101.07197.
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
 *   n_attr_interaction - number of attribute interaction change stats funcs
 *   change_stats_funcs - array of pointers to change statistics functions
 *                        length is n-n_attr-n_dyadic-n_attr_interaction
 *   lambda_values      - array of lambda values for change stats funcs
 *                        same length as change_stats_funcs
 *   attr_change_stats_funcs - array of pointers to change statistics functions
 *                             length is n_attr
 *   dyadic_change_stats_funcs - array of pointers to dyadic change stats funcs
 *                             length is n_dyadic
 *   attr_interaction_change_stats_funcs - array of points to attribute
 *                             interaction change stats funcs. Length is
 *                             n_attr_interaction.
 *   attr_indices   - array of n_attr attribute indices (index into g->binattr
 *                    or g->catattr) corresponding to attr_change_stats_funcs
 *                    E.g. for Sender effect on the first binary attribute,
 *                    attr_indices[x] = 0 and attr_change_stats_funcs[x] =
 *                    changeSender
 *   attr_interaction_pair_indices - array of n_attr_interaction attribute
 *                    pair indices (as above but each element is pair of
 *                    such indices) for attribute interaction effects.
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
 *   forbidReciprocity - if True do not allow reciprocated arcs.
 *   citationERGM      - use cERGM (citation ERGM) estimation conditional
 *                       on term (time period)
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
                    uint_t n_attr_interaction,
                    change_stats_func_t *change_stats_funcs[],
                    double               lambda_values[],
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
                    bool useConditionalEstimation,
                    bool forbidReciprocity, bool citationERGM)
{
  uint_t accepted = 0;    /* number of accepted moves */
  double acceptance_rate;
  uint_t i,j,k,l;
  bool   isDelete = FALSE; /* only init to fix warning */
  double *changestats = (double *)safe_malloc(n*sizeof(double));
  double total;  /* sum of theta*changestats */

  bool allowLoops = TRUE; /* XXX */

  assert(!(citationERGM && useConditionalEstimation)); /* cannot do both */
  assert(!(allowLoops && (useConditionalEstimation || citationERGM))); /* no loops for snowball sampling or citation ERGM */
  

  for (i = 0; i < n; i++)
    addChangeStats[i] = delChangeStats[i] = 0;

  for (k = 0; k < sampler_m; k++) {

    if (useConditionalEstimation) {
      /* Select two nodes i, j in inner waves (i.e. fixing ties in
         outermost wave and between outermost and second-outermost
         waves) uniformly at random, and toggle arc between them,
         however with the additional constraints for snowball sample
         conditional estimation that a tie can only be between adjacent
         waves and a tie cannot be deleted if it is last reamainign tie
         connecting node to preceding wave. Note ignoring arc direction
         here as assumed snowball sample ignored arc directions. */
      assert(!forbidReciprocity); /* TODO not implemented for snowball */
      assert(!allowLoops);
      do {
        i = g->inner_nodes[int_urand(g->num_inner_nodes)];
        do {
          j = g->inner_nodes[int_urand(g->num_inner_nodes)];        
        } while (i == j);
        assert(g->zone[i] < g->max_zone && g->zone[j] < g->max_zone);
        /* any tie must be within same zone or between adjacent zones */
        assert(labs((long)g->zone[i] - (long)g->zone[j]) <= 1 ||
               !isArcIgnoreDirection(g, i, j));
      } while (labs((long)g->zone[i] - (long)g->zone[j]) > 1 ||
               (isArcIgnoreDirection(g, i, j) &&
                ((g->zone[i] > g->zone[j] && g->prev_wave_degree[i] == 1) ||
                 (g->zone[j] > g->zone[i] && g->prev_wave_degree[j] == 1))));
    } else if (citationERGM) {
      /* cERGM: select random node i in last time period (term) and random
       * node j (in any term) and toggle arc i->j between them. In this way
       * we have all arcs (citations) in terms earlier than the last fixed,
       * and we only create citations from nodes in the last term. */
      assert(!allowLoops);
      do {
	i = g->maxterm_nodes[int_urand(g->num_maxterm_nodes)];
	do {
	  j = int_urand(g->num_nodes);
	} while (i == j);
	isDelete = isArc(g, i, j);
      } while (forbidReciprocity && !isDelete && isArc(g, j, i));
      assert(g->term[i] == g->max_term);
    } else {
      /* Basic sampler (no conditional estimation): select two
         nodes i and j uniformly at random and toggle arc between
         them. */
      do {
        i = int_urand(g->num_nodes);
        do {
          j = int_urand(g->num_nodes);
        } while (!allowLoops && i == j);
        isDelete = isArc(g, i ,j);
      } while (forbidReciprocity && !isDelete && isArc(g, j, i));
    }
    
    /* The change statistics are all computed on the basis of adding arc i->j
       so if if the arc exists, we temporarily remove it to compute the
       change statistics, and negate them */
    SAMPLER_DEBUG_PRINT(("%s %d -> %d\n",isDelete ? "del" : "add", i, j));
    if (isDelete) {
      removeArc(g, i, j);
    }

    total = calcChangeStats(g, i, j, n, n_attr, n_dyadic,
                            n_attr_interaction, change_stats_funcs,
                            lambda_values,
                            attr_change_stats_funcs, dyadic_change_stats_funcs,
                            attr_interaction_change_stats_funcs,
                            attr_indices, attr_interaction_pair_indices,
                            theta, isDelete, changestats);
    
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
        insertArc(g, i, j);
      }
    }
  }
  
  acceptance_rate = (double)accepted / sampler_m;
  free(changestats);
  return acceptance_rate;
}
