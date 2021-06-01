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
#include "tntSampler.h"
#include "simulation.h"
#include "loadDigraph.h"

/*****************************************************************************
 *
 * local functions
 *
 ****************************************************************************/


/*
 * Make an Erdos-Renyi aka Bernoulli directed random graph G(n, m)
 * with n nodes and m m edges. I.e. the m edges are placed between
 * dyads chosen uniformly at random.
 *
 * Parameters:
 *   g                        - digraph object. Modifed.
 *   numArcs                  - number of arcs [m in G(n,m) notation]
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
 *   useConditionalEstimation - if True preserve snowball sampling zone 
 *                              structure.
 *   forbidReciprocity        - if True do not allow reciprocated arcs.
 *   addChangeStats           - (Out) vector of n change stats for add moves
 *                              Allocated by caller. Assuming we start
 *                              with the empty graph, this will then
 *                              have the counts of statistics in the 
 *                              graph at the end.
 *  theta                     - parameter values (required by 
 *                              calcChangeStats for total but value
 *                              not used here)
 *
 * Return value:
 *   None. The digraph parameter g is updated.
 */
static void make_erdos_renyi_digraph(digraph_t *g, uint_t numArcs,
                                     uint_t n, uint_t n_attr, uint_t n_dyadic,
                                     uint_t n_attr_interaction,
                                     change_stats_func_t *change_stats_funcs[],
                                     double lambda_values[],
                                     attr_change_stats_func_t
                                                    *attr_change_stats_funcs[],
                                     dyadic_change_stats_func_t
                                                 *dyadic_change_stats_funcs[],
                                     attr_interaction_change_stats_func_t
                                     *attr_interaction_change_stats_funcs[],
                                     uint_t attr_indices[],
                                     uint_pair_t attr_interaction_pair_indices[],                                     
                                     bool useConditionalEstimation,
                                     bool forbidReciprocity,
                                     double addChangeStats[], double theta[])
{
  uint_t i, j, k, l;
  double *changestats = (double *)safe_malloc(n*sizeof(double));
  
  for (k = 0; k < numArcs; k++) {
    if (useConditionalEstimation) {
      assert(!forbidReciprocity); /* TODO not implemented for snowball */
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
    } else {
      /* not using conditional estimation */
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

    /* add change statistics to addChangeStats array */
    (void)calcChangeStats(g, i, j, n, n_attr, n_dyadic, n_attr_interaction,
                          change_stats_funcs,
                          lambda_values,
                          attr_change_stats_funcs,
                          dyadic_change_stats_funcs,
                          attr_interaction_change_stats_funcs,
                          attr_indices,
                          attr_interaction_pair_indices,
                          theta, FALSE, changestats);
    for (l = 0; l < n; l++)
      addChangeStats[l] += changestats[l];
    
    /* actually do the move by adding the arc */
    if (useConditionalEstimation) {
      insertArc_allinnerarcs(g, i, j);
    } else {
      insertArc_allarcs(g, i, j);
    }
  }
  free(changestats);
}


/*****************************************************************************
 *
 * external functions
 *
 ****************************************************************************/


/*
 * Generate digraphs from ERGM distribution with supplied parameters.
 *
 * Parameters:
 *   g      - (in/out) Initial digraph object (empty graph with N nodes
 *            intiially where N is number of nodes in graphs to simulate).
 *   n      - number of parameters (length of theta vector and total
 *            number of change statistics functions)
 *   n_attr - number of attribute change statistics functions
 *   n_dyadic -number of dyadic covariate change stats funcs
 *   n_attr_interaction - number of attribute interaction change stats funcs
 *   change_stats_funcs - array of pointers to change statistics functions
 *                        length is n - n_attr - n_dyadic - n_attr_interaction
 *   lambda_values      - array of lambda values for change stats funcs
 *                        same length as change_stats_funcs
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
 *   sample_size    - number of samples to take from ERGM simulation
 *   interval       - sampler iterations between each sample
 *   burnin         - number of iterations to discard initially 
 *   theta          - array of n parameter values corresponding to
 *                    change stats funcs. Allocated by caller.
 *                    iteration, not just every outer iteration.
 *   useIFDsampler  - if true, use the IFD sampler instead of the basic 
 *                    sampler
 *   ifd_K          - consant for multiplying IFD auxiliary parameter
 *                    (only used if useIFDsampler is True).
 *   useConditionalSimulation - if True, do conditional simulation of 
 *                              snowball network samples.
 *   forbidReciprocity - if True, constrain ERGM sampling so that reciprocated
 *                       arcs are not allowed to be created (so simulation
 *                       is conditional on no reciprocated arcs, should have
 *                       none in input observed graph).
 *   sim_net_file_prefix -  simulated network output filename prefix 
 *   dzA_outfile         - open (write) file to write dzA values to.
 *   outputSimulatedNetworks - if True write simulated networks in Pajek format.
 *   arc_param_index     - index in theta[] parameter of Arc parameter value.
 *                         Only used for useIFDsampler=TRUE
 *   dzA               - (in/Out) vector of n change stats
 *                             Allocated by caller, set to initial graph values
 *   useTNTsampler     - use TNT sampler not IFD or basic.
 *   citationERGM      - use cERGM (citation ERGM) estimation conditional
 *                       on term (time period)
 *
 * Return value:
 *   Nonzero on error, 0 if OK.

 *
 */
int simulate_ergm(digraph_t *g, uint_t n, uint_t n_attr, uint_t n_dyadic,
                  uint_t n_attr_interaction,
                  change_stats_func_t *change_stats_funcs[],
                  double lambda_values[],
                  attr_change_stats_func_t *attr_change_stats_funcs[],
                  dyadic_change_stats_func_t *dyadic_change_stats_funcs[],
                  attr_interaction_change_stats_func_t
                  *attr_interaction_change_stats_funcs[],
                  uint_t attr_indices[],
                  uint_pair_t attr_interaction_pair_indices[],
                  uint_t sample_size, uint_t interval, uint_t burnin,
                  double theta[],
                  bool useIFDsampler, double ifd_K,
                  bool useConditionalSimulation,
                  bool forbidReciprocity,
                  char *sim_net_file_prefix,
                  FILE *dzA_outfile,
                  bool outputSimulatedNetworks,
                  uint_t arc_param_index,
                  double dzA[],
		  bool useTNTsampler,
		  bool citationERGM)
{
  FILE          *sim_outfile;
  char           sim_outfilename[PATH_MAX+1];
  double acceptance_rate = 0;
  double *addChangeStats = (double *)safe_malloc(n*sizeof(double));
  double *delChangeStats = (double *)safe_malloc(n*sizeof(double));
  double dzArc; /* only used for IFD sampler */
  double ifd_aux_param;  /* auxiliary parameter for IFD sampler */
  uint_t l;
  uint_t      samplenum;
  ulonglong_t iternum;
  struct timeval start_timeval, end_timeval, elapsed_timeval;
  int            etime;
  char           suffix[16]; /* only has to be large enough for "_x.txt" 
                                where fx is iteration number */

  assert(!(useIFDsampler && useTNTsampler));
  assert(!(citationERGM && useConditionalSimulation));  
  
  if (useIFDsampler)
    ifd_aux_param = theta[arc_param_index] + arcCorrection(g);
  

  printf("sampleSize = %u, interval = %u burnin = %u\n",
         sample_size, interval, burnin);
  if (useIFDsampler)
    printf("IFD sampler ifd_K = %g initial auxiliary parameter V = %g\n",
           ifd_K, ifd_aux_param);
  else if (useTNTsampler)
    printf("TNT sampler\n");
  if (useConditionalSimulation)
    printf("Doing conditional simulation of snowball sample\n");
  if (forbidReciprocity)
    printf("Simulation is conditional on no reciprocated arcs\n");
  if (citationERGM)
    printf("citation ERGM (cERGM) simulation conditional on term\n");

  if (burnin > 0) {
    gettimeofday(&start_timeval, NULL);
    if (useIFDsampler) {
      assert(!citationERGM); /* TODO implement cERGM in IFD sampler */
      acceptance_rate = ifdSampler(g, n, n_attr, n_dyadic, n_attr_interaction,
                                   change_stats_funcs,
                                   lambda_values,
                                   attr_change_stats_funcs,
                                   dyadic_change_stats_funcs,
                                   attr_interaction_change_stats_funcs,
                                   attr_indices,
                                   attr_interaction_pair_indices,
                                   theta,
                                   addChangeStats, delChangeStats, burnin,
                                   TRUE, /*actually do moves */
                                   ifd_K, &dzArc, &ifd_aux_param,
                                   useConditionalSimulation,
                                   forbidReciprocity);
    } else if (useTNTsampler) {
      assert(!citationERGM); /* TODO implement cERGM in TNT sampler */
      acceptance_rate = tntSampler(g, n, n_attr, n_dyadic,
				   n_attr_interaction,
				   change_stats_funcs,
                                   lambda_values,
				   attr_change_stats_funcs,
				   dyadic_change_stats_funcs,
				   attr_interaction_change_stats_funcs,
				   attr_indices,
				   attr_interaction_pair_indices,
				   theta,
				   addChangeStats, delChangeStats,
				   burnin,
				   TRUE,/*actually do moves*/
				   useConditionalSimulation,
				   forbidReciprocity);
    } else {
      acceptance_rate = basicSampler(g, n, n_attr, n_dyadic,
                                     n_attr_interaction,
                                     change_stats_funcs,
                                     lambda_values,
                                     attr_change_stats_funcs,
                                     dyadic_change_stats_funcs,
                                     attr_interaction_change_stats_funcs,
                                     attr_indices,
                                     attr_interaction_pair_indices,
                                     theta,
                                     addChangeStats, delChangeStats,
                                     burnin,
                                     TRUE,/*actually do moves*/
                                     useConditionalSimulation,
                                     forbidReciprocity,
				     citationERGM);
    }
    for (l = 0; l < n; l++) {
      dzA[l] += addChangeStats[l] - delChangeStats[l]; /* dzA accumulates */
      /* but during burn-in we do not output these values */
    }
    gettimeofday(&end_timeval, NULL);
    timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
    etime = 1000 * elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
    printf("burnin %u iterations took %.2f s\n", burnin, (double)etime/1000);
  }

  for (samplenum = 0; samplenum < sample_size; samplenum++) {
    if (useIFDsampler) {
      assert(!citationERGM); /* TODO implement cERGM in IFD sampler */
      acceptance_rate = ifdSampler(g, n, n_attr, n_dyadic, n_attr_interaction,
                                   change_stats_funcs,
                                   lambda_values,
                                   attr_change_stats_funcs,
                                   dyadic_change_stats_funcs,
                                   attr_interaction_change_stats_funcs,
                                   attr_indices,
                                   attr_interaction_pair_indices,
                                   theta,
                                   addChangeStats, delChangeStats, interval,
                                   TRUE, /*actually do moves */
                                   ifd_K, &dzArc, &ifd_aux_param,
                                   useConditionalSimulation,
                                   forbidReciprocity);
    } else if (useTNTsampler) {
      assert(!citationERGM); /* TODO implement cERGM in TNT sampler */
      acceptance_rate = tntSampler(g, n, n_attr, n_dyadic,
				   n_attr_interaction,
				   change_stats_funcs,
                                   lambda_values,
				   attr_change_stats_funcs,
				   dyadic_change_stats_funcs,
				   attr_interaction_change_stats_funcs,
				   attr_indices,
				   attr_interaction_pair_indices,
				   theta,
				   addChangeStats, delChangeStats,
                                     interval,
				   TRUE,/*actually do moves*/
				   useConditionalSimulation,
				   forbidReciprocity);
    } else {
      acceptance_rate = basicSampler(g, n, n_attr, n_dyadic,
                                     n_attr_interaction,
                                     change_stats_funcs,
                                     lambda_values,
                                     attr_change_stats_funcs,
                                     dyadic_change_stats_funcs,
                                     attr_interaction_change_stats_funcs,
                                     attr_indices,
                                     attr_interaction_pair_indices,
                                     theta,
                                     addChangeStats, delChangeStats,
                                     interval,
                                     TRUE,/*actually do moves*/
                                     useConditionalSimulation,
                                     forbidReciprocity,
				     citationERGM);
    }
    iternum = burnin + interval*(samplenum+1);
    fprintf(dzA_outfile, "%llu ", iternum);
    for (l = 0; l < n; l++) {
      dzA[l] += addChangeStats[l] - delChangeStats[l]; /* dzA accumulates */
      fprintf(dzA_outfile, "%g ", dzA[l]);
    }
    fprintf(dzA_outfile, "%g\n", acceptance_rate);
    fflush(dzA_outfile);

    if (outputSimulatedNetworks) {
      strncpy(sim_outfilename, sim_net_file_prefix,
              sizeof(sim_outfilename)-1);
      sprintf(suffix, "_%llu.net", iternum);
      strncat(sim_outfilename, suffix, sizeof(sim_outfilename) - 1 -
              strlen(suffix));
      sim_outfile = fopen(sim_outfilename, "w");
      write_digraph_arclist_to_file(sim_outfile, g);
      fclose(sim_outfile);
    }
  }
  
  fprintf(stdout, "acceptance rate = %g\n", acceptance_rate);

  free(addChangeStats);
  free(delChangeStats);
  
  return 0;
}

/*
 * Do simulation process using basic or IFD or TNT sampler to draw samples
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
  digraph_t     *g;
  uint_t         n_struct, n_attr, n_dyadic, n_attr_interaction, num_param;
  double        *theta;
  uint_t         i, theta_i, j, k, l, m;
  FILE          *dzA_outfile;
#define HEADER_MAX 65536
  char           fileheader[HEADER_MAX];
  bool           foundArc        = FALSE;
  uint_t         arc_param_index = 0;
  double        *dzA             = NULL;
  FILE          *arclist_file;
  uint_t         num_nodes       = 0;
  double        *changeStats     = NULL;
  
    

  if (!config->stats_filename) {
    fprintf(stderr, "ERROR: statistics output filename statsFile not set\n");
    return -1;
  }
  
  g = allocate_digraph(config->numNodes);
  if (load_attributes(g, config->binattr_filename,
                      config->catattr_filename,
                      config->contattr_filename,
                      config->setattr_filename)) {
    fprintf(stderr, "ERROR: loading node attributes failed\n");
    return -1;
  }

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
  if (config->term_filename) {
    if (add_cergm_terms_to_digraph(g, config->term_filename)) {
      fprintf(stderr, "ERROR: reading cERGM terms from %s failed\n",
              config->term_filename);
      return -1;
    }
#ifdef DEBUG_CERGM
    dump_term_info(g);
#endif /* DEBUG_CERGM */
  }
  
  /* now that we have attributes loaded in g, build the attr_indices
     array in the config struct */
  if (build_attr_indices_from_names(&config->param_config, g) != 0)  {
    fprintf(stderr, "ERROR in attribute parameters\n");
    return -1;
  }
  /* and similary for dyadic covariates */
  if (build_dyadic_indices_from_names(&config->param_config, g, TRUE) != 0)  {
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
    
   /* Ensure that if conditional simulation is to be used, the snowball
      sampling zone structure was specified */
   if (config->useConditionalSimulation) {
     if (!config->zone_filename) {
       fprintf(stderr,
           "ERROR: conditional simulation requested but no zones specified\n");
       return -1;
     }
     if (g->max_zone < 1) {
       fprintf(stderr,
               "ERROR: conditional simulation requested but only one zone\n");
       return -1;
     }
   }  else {
     if (config->zone_filename)
       fprintf(stderr, "WARNING: snowball sampling zones are specified"
	       " but conditional simuilation is not being used\n");
   }


   /* Ensure that if citation ERGM cERGM is to be used, the time period (term)
      values were specified, and also a digraph to read for the initial
      network (from which arcs sent from nodes not in the last term (time
      period) are fixed */
   if (config->citationERGM) {
     if (config->useConditionalSimulation) {
       fprintf(stderr, "ERROR: cannot use both snowball sample conditional"
	       " simulation and citation ERGM\n");
	 return -1;
     }
     if (!config->term_filename) {
       fprintf(stderr,
           "ERROR: citation ERGM simulation requested but no term file\n");
       return -1;
     }
     if (g->max_term < 1) {
       fprintf(stderr,
               "ERROR: citation ERGM simulation requested but only one time period\n");
       return -1;
     }
   }
   
   
   /* 
    *set parameter values from the configuration settings and write
    *  parameters and their values to stdout 
    */
   theta = (double *)safe_calloc(num_param, sizeof(double));
   theta_i = 0;
   printf("\n");
   for (i = 0; i < config->param_config.num_change_stats_funcs; i++, theta_i++){
     theta[theta_i] = config->param_config.param_values[i];
     if (config->param_config.param_lambdas[i] > 0.0) {
       printf("%s(%g) = %g\n", config->param_config.param_names[i],
              config->param_config.param_lambdas[i],
              theta[theta_i]);
     } else {
       printf("%s = %g\n", config->param_config.param_names[i], theta[theta_i]);     }
   }
   
   for (i = 0; i < config->param_config.num_attr_change_stats_funcs;
        i++, theta_i++)  {
     theta[theta_i] = config->param_config.attr_param_values[i];
     printf("%s_%s = %g\n", config->param_config.attr_param_names[i],
            config->param_config.attr_names[i], theta[theta_i]);
   }
   
   for (i = 0; i < config->param_config.num_dyadic_change_stats_funcs;
        i++, theta_i++) {
     theta[theta_i] = config->param_config.dyadic_param_values[i];
     printf("%s = %g\n", config->param_config.dyadic_param_names[i],
            theta[theta_i]);
   }
   
   for (i = 0; i < config->param_config.num_attr_interaction_change_stats_funcs;
        i++, theta_i++)  {
     fprintf(stderr, "TODO: initial values not implemented for interaction effects yet\n");assert(FALSE);     
     /*TODO: theta_theta[i] = config->param_config.attr_interaction_param_values[i]; */
     printf("%s_%s_%s = %g\n",
            config->param_config.attr_interaction_param_names[i],
            config->param_config.attr_interaction_pair_names[i].first,
            config->param_config.attr_interaction_pair_names[i].second,
            theta[theta_i]);
   }
   printf("\n");

   /* Only one sampler can be used (only binary attributes in config,
      did not include multiple options (maybe should) */
   if (config->useIFDsampler && config->useTNTsampler) {
     fprintf(stderr, "ERROR: Only one of the useIFDsampler and"
	     " useTNTsampler options may be used\n");
     return -1;
   }

   /* Ensure that for the IFD sampler the  Arc parameter included, and
      get its index as it is needed to compute the initial value of the
      IFD sampler auxilliary parameter,
      and also instead the numArcs parameter is included (fixed density) */
   if (config->useIFDsampler) {
     for (i = 0; i < config->param_config.num_change_stats_funcs; i++) {
       if (strcasecmp(config->param_config.param_names[i], ARC_PARAM_STR) == 0) {
         arc_param_index = i;
         foundArc = TRUE;
       }
     }
     if (!foundArc) {
       fprintf(stderr, 
               "ERROR: must include %s parameter when using IFD sampler.\n",
               ARC_PARAM_STR);
       return -1;
     }
     if (config->numArcs == 0) {
       fprintf(stderr, "ERROR: must specify nonzero numArcs when "
               "using IFD sampler\n");
       return -1;
     }
   } else  {
     if (config->numArcs != 0) {
       fprintf(stderr, "ERROR: numArcs only used for IFD sampler\n");
       return -1;
     }
   }

   /* Give warnings if parameters set that are not used in selected
      algorithm variation */
   if (!config->useIFDsampler &&
       !DOUBLE_APPROX_EQ(config->ifd_K, DEFAULT_IFD_K)) {
     fprintf(stderr,
             "WARNING: ifd_K is set to %g not default value"
             " but IFD sampler not used\n", config->ifd_K);
   }


   /* allocate change statistics array  */
   dzA = (double *)safe_calloc(num_param, sizeof(double));
   /* set values of graph stats for empty graph; most (but not all) are zero */
   empty_graph_stats(g, num_param, n_attr, n_dyadic,
		     n_attr_interaction,
		     config->param_config.change_stats_funcs,
		     config->param_config.param_lambdas,
		     config->param_config.attr_change_stats_funcs,
		     config->param_config.dyadic_change_stats_funcs,
		     config->param_config.attr_interaction_change_stats_funcs,
		     config->param_config.attr_indices,
		     config->param_config.attr_interaction_pair_indices,
		     dzA);
   if (config->citationERGM) {
     /* For citation ERGM, initialize the graph from the observed
	graph specified in the Pajek format arclist file. All arcs
	sent from nodes in time periods other than the last will be
	fixed; arcs from the last time period will be deleted in the
	initialization here, and can be added and deleted in the
	simulation */
     if (!config->arclist_filename) {
       fprintf(stderr, "ERROR: ctiation ERGM simulation requested but no arclistFile specified.\n");
       return -1;
     }
     if (!(arclist_file = fopen(config->arclist_filename, "r"))) {
       fprintf(stderr, "error opening file %s (%s)\n", 
	       config->arclist_filename, strerror(errno));
       return -1;
     }
     if (config->numArcs != 0) {
       fprintf(stderr, "WARNING: numArcs is set to %u but using citationERGM"
	       " so numArcs parameter is ignored\n", config->numArcs);
     }
     num_nodes = get_num_vertices_from_arclist_file(arclist_file);/* closes file */
     if (num_nodes != config->numNodes) {
       fprintf(stderr, "ERROR: number of nodes specified in config file (%u) does not match number in Pajek file %s (%u)\n", config->numNodes, config->arclist_filename, num_nodes);
       return 1;
     }
     if (!(arclist_file = fopen(config->arclist_filename, "r"))) {
       fprintf(stderr, "error opening file %s (%s)\n", 
	       config->arclist_filename, strerror(errno));
       return -1;
     }
     g = load_digraph_from_arclist_file(arclist_file, g,
					TRUE /*computeStats*/,
					num_param,
					n_attr, n_dyadic, n_attr_interaction,
					config->param_config.change_stats_funcs,
					config->param_config.param_lambdas,
					config->param_config.attr_change_stats_funcs,
					config->param_config.dyadic_change_stats_funcs,
					config->param_config.attr_interaction_change_stats_funcs,
					config->param_config.attr_indices,
					config->param_config.attr_interaction_pair_indices,
					dzA, theta);
     /* Now delete all arcs sent from nodes in the last time period
	(term) so that we start the simulation with an empty set of
	these non-fixed potential arcs (all the other arcs and
	non-arcs, from nodes in previous terms, are fixed). */
     changeStats = (double *)safe_calloc(num_param, sizeof(double));     
     for (k = 0; k < g->num_maxterm_nodes; k++) {
       i = g->maxterm_nodes[k];
       assert(g->term[i] == g->max_term);
       for (l = 0; l < g->outdegree[i]; l++) {
	 j = g->arclist[i][l];
	 assert(isArc(g, i, j));
	 removeArc(g, i, j);
	 /* Update the statistics for removing this arc */
	 calcChangeStats(g, i, j,
			 num_param, n_attr, n_dyadic,
			 n_attr_interaction,
			 config->param_config.change_stats_funcs,
			 config->param_config.param_lambdas,
			 config->param_config.attr_change_stats_funcs,
			 config->param_config.dyadic_change_stats_funcs,
			 config->param_config.attr_interaction_change_stats_funcs,
			 config->param_config.attr_indices,
			 config->param_config.attr_interaction_pair_indices,
			 theta,
			 TRUE, /*isDelete*/
			 changeStats);
	 for (m = 0; m < num_param; m++) {
	   dzA[m] += changeStats[m]; /* isDelete=TRUE above makes change -ve */
	 }
       }
     }
   } else {
     if (config->useIFDsampler) {
       /* Initialize the graph to random (E-R aka Bernoulli) graph with
	  specified number of arcs for fixed density simulation (IFD sampler),
	  and also for TNT sampler since it does 50% add/delete moves */
       make_erdos_renyi_digraph(g, config->numArcs,
				num_param, n_attr, n_dyadic, n_attr_interaction,
				config->param_config.change_stats_funcs,
				config->param_config.param_lambdas,
				config->param_config.attr_change_stats_funcs,
				config->param_config.dyadic_change_stats_funcs,
				config->param_config.attr_interaction_change_stats_funcs,
				config->param_config.attr_indices,
				config->param_config.attr_interaction_pair_indices,                              
				config->useConditionalSimulation,
				config->forbidReciprocity,
				dzA, theta);
     } else if (config->numArcs != 0) {
       fprintf(stderr, "WARNING: numArcs is set to %u but not using IFD sampler"
	       " so numArcs parameter is ignored\n", config->numArcs);
     }
   }

   /* Open the output file for writing */
   if (!(dzA_outfile = fopen(config->stats_filename, "w"))) {
     fprintf(stderr, "ERROR: could not open file %s for writing "
             "(%s)\n", config->stats_filename, strerror(errno));
     return -1;
   }

   /* write headers for statistics output file */
  sprintf(fileheader, "t");
  /* Print the lambda (decay) [hyper-]parameter value for parameters which
     use it (i.e. for the "alternating" statistics); it is 0 for those
     for which it is not applicable.
     Format is to put it in parens after the name e.g. AltTwoPathsTD(2.5) */
  for (i = 0; i < config->param_config.num_change_stats_funcs; i++) {
    if (config->param_config.param_lambdas[i] > 0.0) {
      snprintf(fileheader+strlen(fileheader), HEADER_MAX," %s(%g)",
               config->param_config.param_names[i],
               config->param_config.param_lambdas[i]);
    } else {
      snprintf(fileheader+strlen(fileheader), HEADER_MAX," %s",
               config->param_config.param_names[i]);      
    }
  }
  
  for (i = 0; i < config->param_config.num_attr_change_stats_funcs; i++) 
    snprintf(fileheader+strlen(fileheader), HEADER_MAX, " %s_%s",
             config->param_config.attr_param_names[i],
             config->param_config.attr_names[i]);
  
   for (i = 0; i < config->param_config.num_dyadic_change_stats_funcs; i++)
     snprintf(fileheader+strlen(fileheader), HEADER_MAX, " %s",
              config->param_config.dyadic_param_names[i]);

   for (i = 0; i < config->param_config.num_attr_interaction_change_stats_funcs; i++) 
     snprintf(fileheader+strlen(fileheader), HEADER_MAX, " %s_%s_%s",
              config->param_config.attr_interaction_param_names[i],
              config->param_config.attr_interaction_pair_names[i].first,
              config->param_config.attr_interaction_pair_names[i].second);

   fprintf(dzA_outfile,  "%s AcceptanceRate\n", fileheader);

   print_data_summary(g);
   print_zone_summary(g);
   print_term_summary(g);

#ifdef DEBUG_SIMULATE
   printf("initial graph stats: ");   
   for(i = 0; i < num_param; i++) {
     printf("%g ", dzA[i]);
   }
   printf("\n");
#endif


   printf("\nrunning simulation...\n");
   gettimeofday(&start_timeval, NULL);
   
   simulate_ergm(g, num_param, n_attr, n_dyadic, n_attr_interaction,
                 config->param_config.change_stats_funcs,
                 config->param_config.param_lambdas,
                 config->param_config.attr_change_stats_funcs,
                 config->param_config.dyadic_change_stats_funcs,
                 config->param_config.attr_interaction_change_stats_funcs,
                 config->param_config.attr_indices,
                 config->param_config.attr_interaction_pair_indices,
                 config->sampleSize, config->interval, config->burnin,
                 theta,
                 config->useIFDsampler, config->ifd_K,
                 config->useConditionalSimulation,
                 config->forbidReciprocity,
                 config->sim_net_file_prefix,
                 dzA_outfile,
                 config->outputSimulatedNetworks, arc_param_index,
                 dzA, config->useTNTsampler, config->citationERGM);

   gettimeofday(&end_timeval, NULL);
   timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
   etime = 1000 * elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
   printf("simulation took %.2f s\n", (double)etime/1000);

   fclose(dzA_outfile);

   print_data_summary(g);
     
   free(theta);
   free(dzA);
   free_digraph(g);
   
  return 0;
}
