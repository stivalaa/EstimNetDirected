#ifndef SIMULATION_H
#define SIMULATION_H
/*****************************************************************************
 * File:    simulation.c
 * Author:  Alex Stivala
 * Created: October 2019
 *
 * Draw samples from ERGM distribution of directed graphs.
 *
 ****************************************************************************/

#include "simconfigparser.h"
#include "changeStatisticsDirected.h"

int do_simulation(sim_config_t *config);


#endif /* SIMULATION_H */

