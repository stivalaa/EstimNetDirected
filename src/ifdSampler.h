#ifndef IFDSAMPLER_H
#define IFDSAMPLER_H
/*****************************************************************************
 * 
 * File:    ifdSampler.h
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
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

#include "changeStatisticsDirected.h"

double arcCorrection(const digraph_t *g);

double ifdSampler(digraph_t *g,  uint_t n, uint_t n_attr, uint_t n_dyadic,
                  change_stats_func_t *change_stats_funcs[],
                  attr_change_stats_func_t *attr_change_stats_funcs[],
                  dyadic_change_stats_func_t *dyadic_change_stats_funcs[],
                  uint_t attr_indices[], double theta[],
                  double addChangeStats[], double delChangeStats[], 
                  uint_t sampler_m,
                  bool performMove,
                  double ifd_K, 
                  double *dzArc, double *ifd_aux_param);

#endif /* IFDSAMPLER_H */

