#ifndef SIMULATION_H
#define SIMULATION_H
/*****************************************************************************
 * File:    simulation.c
 * Author:  Alex Stivala
 * Created: October 2019
 *
 * Draw samples from ERGM distribution of graphs.
 *
 ****************************************************************************/

#include "simconfigparser.h"
#include "changeStatisticsDirected.h"

int simulate_ergm(graph_t *g, uint_t n, uint_t n_attr, uint_t n_dyadic,
                  uint_t n_attr_interaction,
                  change_stats_func_t *change_stats_funcs[],
                  double lambda_values[],
                  attr_change_stats_func_t *attr_change_stats_funcs[],
                  dyadic_change_stats_func_t *dyadic_change_stats_funcs[],
                  attr_interaction_change_stats_func_t
                  *attr_interaction_change_stats_funcs[],
                  uint_t attr_indices[],
                  double exponent_values[],
                  uint_pair_t attr_interaction_pair_indices[],
                  uint_t sample_size, ulong_t interval, ulong_t burnin,
                  double theta[],
                  bool useIFDsampler, double ifd_K,
                  bool useConditionalEstimation,
                  bool forbidReciprocity,
                  char *sim_net_file_prefix,
                  FILE *dzA_outfile,
                  bool outputSimulatedNetworks,
                  uint_t arc_param_index,
                  double addChangeStats[], bool useTNTsampler,
                  bool citationERGM, bool allowLoops);

int do_simulation(sim_config_t *config);


#endif /* SIMULATION_H */

