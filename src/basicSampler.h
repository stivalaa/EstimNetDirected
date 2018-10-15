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
 *
 ****************************************************************************/

#include "changeStatisticsDirected.h"

double basicSampler(digraph_t *g,  uint_t n, uint_t n_attr, uint_t n_dyadic,
                    change_stats_func_t *change_stats_funcs[],
                    attr_change_stats_func_t *attr_change_stats_funcs[],
                    dyadic_change_stats_func_t *dyadic_change_stats_funcs[],
                    uint_t attr_indices[], double theta[],
                    double addChangeStats[], double delChangeStats[], 
                    uint_t sampler_m,
                    bool performMove,
                    bool useConditionalEstimation);

#endif /* BASICSAMPLER_H */

