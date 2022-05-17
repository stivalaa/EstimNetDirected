#ifndef BASICSAMPLER_H
#define BASICSAMPLER_H
/*****************************************************************************
 * 
 * File:    basicSampler.h
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

#include "changeStatisticsDirected.h"

double basicSampler(graph_t *g,  uint_t n, uint_t n_attr, uint_t n_dyadic,
                    uint_t n_attr_interaction,
                    change_stats_func_t *change_stats_funcs[],
                    double lambda_values[],
                    attr_change_stats_func_t *attr_change_stats_funcs[],
                    dyadic_change_stats_func_t *dyadic_change_stats_funcs[],
                    attr_interaction_change_stats_func_t
                                     *attr_interaction_change_stats_funcs[],
                    uint_t attr_indices[],
                    uint_pair_t attr_interaction_pair_indices[],
                    double theta[],
                    double addChangeStats[], double delChangeStats[],
                    ulong_t sampler_m,
                    bool performMove,
                    bool useConditionalEstimation,
                    bool forbidReciprocity, bool citationERGM,
                    bool allowLoops);

#endif /* BASICSAMPLER_H */

