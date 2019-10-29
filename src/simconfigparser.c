/*****************************************************************************
 * 
 * File:    simconfigparser.c
 * Author:  Alex Stivala
 * Created: October 2019
 *
 * Parse the simulation configuration file to get algorithm
 * paraemeters, input filenames, parameters to estimate, etc.
 *
 * The config file is a text file with comments marked by '#'
 * character, and "keyword = value" pairs, with structural and attribute
 * parameters specified as e.g. "structParams =  {Arc = -4.0, Reciprocity=2.1}",
 * "attrParams = {Sender(binattrname) = 1.2,Receiver(binattrname)=0.5}", etc.
 * See simconfig.txt for example config file.
 *
 * Tokens and attribute names are not case sensitive, only filenames are.
 *
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <stddef.h>
#include <assert.h>
#include "simconfigparser.h"

/*****************************************************************************
 *
 * constant definitions
 *
 ****************************************************************************/


/* 
 * Configuration parameter names, types and descriptions. These 
 * define what is parsed from the config file, and the help/error messages.
 * Configuration parameter names (keywords) are not case sensitive.
 */
const config_param_t SIM_CONFIG_PARAMS[] = {
  {"numNodes",      PARAM_TYPE_UINT,     offsetof(sim_config_t, numNodes),
   "number of nodes in digraph"},
  
  {"sampleSize",  PARAM_TYPE_UINT,     offsetof(sim_config_t, sampleSize),
   "number of network samples to take from simulation"},

  {"interval",  PARAM_TYPE_UINT,     offsetof(sim_config_t, interval),
   "interval (iterations) between samples"},

  {"burnin",  PARAM_TYPE_UINT,     offsetof(sim_config_t, burnin),
   "number of iterations to throw away before first sample"},

  {"useIFDsampler", PARAM_TYPE_BOOL,    offsetof(sim_config_t, useIFDsampler),
   "use Improved Fixed Density sampler instead of basic sampler"},

  {"ifd_K",         PARAM_TYPE_DOUBLE,  offsetof(sim_config_t, ifd_K),
   "multiplier for auxiliary parameter step size in IFD sampler"},

  {"outputSimulatedNetwork", PARAM_TYPE_BOOL,
   offsetof(sim_config_t, outputSimulatedNetwork),
   "output simulated network in Pajek format at end of MCMC simulation"},

  {"binattrFile",   PARAM_TYPE_STRING,   offsetof(sim_config_t, binattr_filename),
  "binary attributes file"},

  {"catattrFile",   PARAM_TYPE_STRING,   offsetof(sim_config_t, catattr_filename),
  "categorical attributes file"},

  {"contattrFile",  PARAM_TYPE_STRING,   offsetof(sim_config_t, contattr_filename),
  "continuous attributes file"},

  {"setattrFile",   PARAM_TYPE_STRING,   offsetof(sim_config_t, setattr_filename),
  "set attributes file"},

  {"statsFile",  PARAM_TYPE_STRING,  offsetof(sim_config_t, stats_filename),
   "statistics output filename"},

  {"simNetFilePrefix", PARAM_TYPE_STRING,
   offsetof(sim_config_t, sim_net_file_prefix),
   "simulated network output file prefix"},

  {"zoneFile",      PARAM_TYPE_STRING,  offsetof(sim_config_t, zone_filename),
   "snowball sample zone file"},

  {"useConditionalEstimation", PARAM_TYPE_BOOL,
   offsetof(sim_config_t, useConditionalEstimation),
   "do conditional estimation for snowball network sample"},

  {"forbidReciprocity",PARAM_TYPE_BOOL, offsetof(sim_config_t, forbidReciprocity),
   "constrain ERGM sampler to not allow reciprocated arcs"},
  /* This is useful for graphs that have no reciprocated arcs in the observed
     graph so cannot use Reciprocity parameter, and want to enforce 
     constraint that there can be no reciprocated arcs e.g. for a citatoin
     network. 
     TODO should have some more general way of specifying constraints
     like ergm-constraints in statnet instead of this ad-hoc way */


  {STRUCT_PARAMS_STR,  PARAM_TYPE_SET,      0, /*no offset, coded explicitly*/
  "structural parameters to estimate"},

  {ATTR_PARAMS_STR, PARAM_TYPE_SET,         0, /*no offset, coded explicitly*/
  "binary/categorical/continuous/set attribute parameters to estimate"},

  {DYADIC_PARAMS_STR, PARAM_TYPE_SET,       0, /*no offset, coded explicitly*/
  "dyadic covariate parameters to estimate"},

  {ATTR_INTERACTION_PARAMS_STR,PARAM_TYPE_SET, 0,/*no offset, coded explicitly*/
   "attribute pair interaction parameters to estimate"}
};
const uint_t NUM_SIM_CONFIG_PARAMS = sizeof(SIM_CONFIG_PARAMS) /
  sizeof(SIM_CONFIG_PARAMS[0]);



/*****************************************************************************
 *
 * externally visible variables
 *
 ****************************************************************************/

/* 
 * The config structure here is initialized with default values, and
 * check_and_set_param_value() updates it with values from the config file.
 * String values must be NULL (not string constant) so that they can
 * be free()ed later (if a value is specified). The default values (if any) are 
 * therefore set in code with strdup() just like specified values.
 */
sim_config_t SIM_CONFIG = {
  0,     /* numNodes */
  SIM_DEFAULT_SAMPLE_SIZE,/* sampleSize */
  SIM_DEFAULT_INTERVAL,   /* interval */
  SIM_DEFAULT_BURNIN,     /* burnin */
  FALSE, /* useIFDsampler */
  SIM_DEFAULT_IFD_K,   /* ifd_K */
  FALSE, /* outputSimulatedNetwork */
  NULL,  /* binattr_filename */
  NULL,  /* catattr_filename */
  NULL,  /* contattr_filename */
  NULL,  /* setattr_filename */
  NULL,  /* stats_filename */
  NULL,  /* sim_net_file_prefix */
  NULL,  /* zone_filename */
  FALSE, /* useConditionalEstimation */
  FALSE, /* forbidReciprocity */
  {
    0,     /* num_change_stats_funcs */
    NULL,  /* change_stats_funcs */
    NULL,  /* param_names */
    0,     /* num_attr_change_stats_funcs */
    NULL,  /* attr_change_stats_funcs */
    NULL,  /* attr_names */
    NULL,  /* attr_indices */
    NULL,  /* attr_param_names */
    0,     /* num_dyadic_change_stats_funcs */
    NULL,  /* dyadic_change_stats_funcs */
    NULL,  /* dyadic_names */
    NULL,  /* dyadic_indices */
    NULL,  /* dyadic_types */
    NULL,  /* dyadic_param_names */
    0,     /* num_attr_interaction_change_stats_funcs */
    NULL,  /* attr_attr_interaction_change_stats_funcs */
    NULL,  /* attr_interaction_names */
    NULL,  /* attr_interaction_indices */
    NULL   /* attr_interaction_param_names */
  } /* param_config */
};


/*****************************************************************************
 *
 * file static variables
 *
 ****************************************************************************/


/* 
 * Array of flags to set to TRUE when a parameter value is set. 
 * If we don't check this then the same parameter can be set multiple
 * time (and the last one takes effect, or in the case of structParams
 * and attrParams, they are just added to the list - this could be fine,
 * but better readability of config file if we give error and only allow
 * a single instance of each parameter).
 * Note that since structParams and attrParams are treated specially, not
 * like the simple 'name = value' parameters, we do not use this,
 * but simply check CONFIG.num_change_stats_funcs and 
 * CONFIG.num_attr_change_stats_funcs
 */
static bool SIM_CONFIG_IS_SET[] = {
  FALSE, /* numNodes */ 
  FALSE, /* sampleSize */
  FALSE, /* interval */
  FALSE, /* burnin */
  FALSE, /* useIFDsampler */
  FALSE, /* ifd_K */
  FALSE, /* outputSimulatedNetwork */
  FALSE, /* binattr_filename */
  FALSE, /* catattr_filename */
  FALSE, /* contattr_filename */
  FALSE, /* setattr_filename */
  FALSE, /* stats_filename */
  FALSE, /* sim_net_file_prefix */
  FALSE, /* zone_filename */
  FALSE, /* useConditionalEstimation */
  FALSE, /* forbidReciprocity */
  FALSE, /* (NOT USED) structParams */
  FALSE, /* (NOT USED) attrParams */
  FALSE, /* (NOT USED) dyadicParams */
  FALSE  /* (NOT USED) attrInteractionParams */
};
static const uint_t NUM_SIM_CONFIG_IS_SET = sizeof(SIM_CONFIG_IS_SET)/sizeof(SIM_CONFIG_IS_SET[0]);



/*****************************************************************************
 *
 * local functions
 *
 ****************************************************************************/






/*****************************************************************************
 *
 * external functions
 *
 ****************************************************************************/

/*
 * Parse the configuration file
 *
 * Parameters:
 *   config_filename - filename of the configuration file to read
 *
 * Return value:
 *   Pointer to structure with parsed confiugration values or NULL on error.
 */
sim_config_t *parse_sim_config_file(const char *config_filename)
{
  char        paramname[TOKSIZE];  /* parameter name buffer */
  char        value[TOKSIZE];      /* parameter value buffer */
  FILE       *config_file;
  int         rc;

  if (!(config_file = fopen(config_filename, "r"))) {
    fprintf(stderr, "ERROR: could not open configuration file %s (%s)\n",
            config_filename, strerror(errno));
    return NULL;
  }
  while (!feof(config_file) && !ferror(config_file)) {
    if ((rc = get_paramname_value(config_file, paramname, value)) == 0) {
      if (check_and_set_param_value(paramname, value, config_file,
                                    &SIM_CONFIG, SIM_CONFIG_IS_SET,
                                    &SIM_CONFIG.param_config,
                                    SIM_CONFIG_PARAMS, NUM_SIM_CONFIG_PARAMS) != 0)
        return NULL;
    } else if (rc < 0)
      return NULL;
  }
  if (ferror(config_file)) {
    fprintf(stderr, "ERROR: error reading configuration file %s (%s)\n",
            config_filename, strerror(errno));
    return NULL;
  }
  fclose(config_file);
  return &SIM_CONFIG; /* return pointer to static CONFIG structure */
}



/*
 * Free the config structure returned by parse_config_file()
 *
 * Parameters:
 *     config - pointer to config struct returned by parse_config_file()
 *
 * Return value:
 *     None
 */
void free_sim_config_struct(sim_config_t *config)
{
  /* In fact parse_config_file() returns pointer to static CONFIG struct,
     so just free the pointers inside it */
  if (!config)
    return;
  assert(config == &SIM_CONFIG);
  free(config->binattr_filename);
  free(config->catattr_filename);
  free(config->contattr_filename);
  free(config->setattr_filename);
  free(config->stats_filename);
  free(config->sim_net_file_prefix);
  free(config->zone_filename);
  free_param_config_struct(&config->param_config);
}



/*
 * Initialize the parser. This consists only of setting the default string
 * parameter values (see comments on initialization of static CONFIG struct)
 */
void init_sim_config_parser(void)
{
  assert(NUM_SIM_CONFIG_IS_SET == NUM_SIM_CONFIG_PARAMS);
  SIM_CONFIG.sim_net_file_prefix = safe_strdup("simulation");
}



