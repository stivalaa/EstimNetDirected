/*****************************************************************************
 * 
 * File:    simulation.c
 * Author:  Alex Stivala
 * Created: October 2019
 *
 * Draw samples from ERGM distribution of directed graphs.
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include "utils.h"
#include "digraph.h"
#include "basicSampler.h"
#include "ifdSampler.h"
#include "simulation.h"


/*
 * Do simulation process using basic or IFD sampler to draw samples
 * from ERGM digraph distribution.
 *
 * Parameters:
 *   config - (in/out)configuration settings structure  - this is modified
 *            by calling build_attr_indices_from_names() etc.
 *
 * Return value:
 *    0 if OK else -ve value for error.
 */
int do_simulation(sim_config_t * config)
{
  struct timeval start_timeval, end_timeval, elapsed_timeval;
  int            etime;
  FILE          *arclist_file;
  digraph_t     *g;
  uint_t         i;
  FILE          *sim_outfile;
  char           sim_outfilename[PATH_MAX+1];
  uint_t         n_struct, n_attr, n_dyadic, n_attr_interaction, num_param;
  
  if (!(arclist_file = fopen(config->arclist_filename, "r"))) {
    fprintf(stderr, "error opening file %s (%s)\n", 
            config->arclist_filename, strerror(errno));
    return -1;
  }
  gettimeofday(&start_timeval, NULL);
  printf("loading arc list from %s and building two-path matrices...",
         config->arclist_filename);
  g = load_digraph_from_arclist_file(arclist_file,
                                     config->binattr_filename,
                                     config->catattr_filename,
                                     config->contattr_filename,
                                     config->setattr_filename);
  gettimeofday(&end_timeval, NULL);
  timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
  etime = 1000 * elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
  printf("%.2f s\n", (double)etime/1000);
#ifdef DEBUG_DIGRAPH
  dump_digraph_arclist(g);
#endif /*DEBUG_DIGRAPH*/

  if (config->zone_filename) {
    if (add_snowball_zones_to_digraph(g, config->zone_filename)) {
      fprintf(stderr, "ERROR: reading snowball sampling zones from %s failed\n",
              config->zone_filename);
      return -1;
    }
#ifdef DEBUG_SNOWBALL
    dump_zone_info(g);
#endif /* DEBUG_SNOWBALL */
  }
  
  
  print_data_summary(g);
  print_zone_summary(g);


  /* now that we have attributes loaded in g, build the attr_indices
     array in the config struct */
  if (build_attr_indices_from_names(&config->param_config, g) != 0)  {
    fprintf(stderr, "ERROR in attribute parameters\n");
    return -1;
  }
  /* and similary for dyadic covariates */
  if (build_dyadic_indices_from_names(&config->param_config, g) != 0)  {
    fprintf(stderr, "ERROR in dyadic covariate parameters\n");
    return -1;
  }
  /* and attribute interaction parameters */
  if (build_attr_interaction_pair_indices_from_names(&config->param_config, g) != 0) {
    fprintf(stderr, "ERROR in attribute interaction parameters\n");
    return -1;
  }

  /* note num_param is computed here as build_dyadic_indices_from_names()
     can decrease config->num_dyadic_change_stats_funcs from its 
     initial value */
  n_struct = config->param_config.num_change_stats_funcs;
  n_attr = config->param_config.num_attr_change_stats_funcs;
  n_dyadic = config->param_config.num_dyadic_change_stats_funcs;
  n_attr_interaction = config->param_config.num_attr_interaction_change_stats_funcs;
  num_param =  n_struct + n_attr + n_dyadic + n_attr_interaction;
    
   /* Ensure that if conditional estimation is to be used, the snowball
      sampling zone structure was specified */
   if (config->useConditionalEstimation) {
     if (!config->zone_filename) {
       fprintf(stderr,
           "ERROR: conditional estimation requested but no zones specified\n");
       return -1;
     }
     if (g->max_zone < 1) {
       fprintf(stderr,
               "ERROR: conditional estimation requested but only one zone\n");
       return -1;
     }
   }

  free_digraph(g);

  return 0;
}
