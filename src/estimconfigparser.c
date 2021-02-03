/*****************************************************************************
 * 
 * File:    estimconfigparser.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * Parse the estimation configuration file to get algorithm
 * paraemeters, input filenames, parameters to estimate, etc.
 *
 * The config file is a text file with comments marked by '#'
 * character, and "keyword = value" pairs, with structural and attribute
 * parameters specified as e.g. "structParams =  {Arc ,  Reciprocity}",
 * "attrParams = {Sender(binattrname),Receiver(binattrname)}", etc.
 * See config.txt for example config file.
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
#include "estimconfigparser.h"

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
const config_param_t ESTIM_CONFIG_PARAMS[] = {
  {"ACA_S",         PARAM_TYPE_DOUBLE,   offsetof(estim_config_t, ACA_S),
   "multiplier for step size in Algorithm S"},

  {"ACA_EE",        PARAM_TYPE_DOUBLE,   offsetof(estim_config_t, ACA_EE),
   "multiplier for step size in Algorithm EE"},

  {"compC",         PARAM_TYPE_DOUBLE,   offsetof(estim_config_t, compC),
   "multiplier of sd(theta)/mean(theta) to limit variance"},

  {"samplerSteps",  PARAM_TYPE_UINT,     offsetof(estim_config_t, samplerSteps),
   "sampler iterations (per algorithm step)"},

  {"Ssteps",        PARAM_TYPE_UINT,     offsetof(estim_config_t, Ssteps),
   "steps of Algorithm S"},

  {"EEsteps",       PARAM_TYPE_UINT,     offsetof(estim_config_t, EEsteps),
   "steps of Algorithm EE"},

  {"EEinnerSteps",  PARAM_TYPE_UINT,     offsetof(estim_config_t, EEinnerSteps),
   "inner iterations of Algorithm EE"},

  {"outputAllSteps", PARAM_TYPE_BOOL,    offsetof(estim_config_t, outputAllSteps),
   "output theta and dzA values on every iteration of EE algorithm)"},

  {"useIFDsampler", PARAM_TYPE_BOOL,    offsetof(estim_config_t, useIFDsampler),
   "use Improved Fixed Density sampler instead of basic sampler"},

  {"useTNTsampler", PARAM_TYPE_BOOL,    offsetof(estim_config_t, useTNTsampler),
   "use Tie-No-Tie sampler instead of basic or IFD sampler"},

  {"ifd_K",         PARAM_TYPE_DOUBLE,  offsetof(estim_config_t, ifd_K),
   "multiplier for auxiliary parameter step size in IFD sampler"},

  {"outputSimulatedNetwork", PARAM_TYPE_BOOL,
   offsetof(estim_config_t, outputSimulatedNetwork),
   "output simulated network in Pajek format at end of MCMC simulation"},

  {"arclistFile",   PARAM_TYPE_STRING,   offsetof(estim_config_t, arclist_filename),
  "Network in Pajek arc list format"},

  {"binattrFile",   PARAM_TYPE_STRING,   offsetof(estim_config_t, binattr_filename),
  "binary attributes file"},

  {"catattrFile",   PARAM_TYPE_STRING,   offsetof(estim_config_t, catattr_filename),
  "categorical attributes file"},

  {"contattrFile",  PARAM_TYPE_STRING,   offsetof(estim_config_t, contattr_filename),
  "continuous attributes file"},

  {"setattrFile",   PARAM_TYPE_STRING,   offsetof(estim_config_t, setattr_filename),
  "set attributes file"},

  {"thetaFilePrefix",PARAM_TYPE_STRING,  offsetof(estim_config_t, theta_file_prefix),
   "theta output file prefix"},

  {"dzAFilePrefix",  PARAM_TYPE_STRING,  offsetof(estim_config_t, dzA_file_prefix),
   "dzA output file prefix"},
  
  {"simNetFilePrefix", PARAM_TYPE_STRING,
   offsetof(estim_config_t, sim_net_file_prefix),
   "simulated network output file prefix"},

  {"zoneFile",      PARAM_TYPE_STRING,  offsetof(estim_config_t, zone_filename),
   "snowball sample zone file"},

  {"useConditionalEstimation", PARAM_TYPE_BOOL,
   offsetof(estim_config_t, useConditionalEstimation),
   "do conditional estimation for snowball network sample"},

  {"forbidReciprocity",PARAM_TYPE_BOOL, offsetof(estim_config_t, forbidReciprocity),
   "constrain ERGM sampler to not allow reciprocated arcs"},
  /* This is useful for graphs that have no reciprocated arcs in the observed
     graph so cannot use Reciprocity parameter, and want to enforce 
     constraint that there can be no reciprocated arcs e.g. for a citatoin
     network. 
     TODO should have some more general way of specifying constraints
     like ergm-constraints in statnet instead of this ad-hoc way */

  {"useBorisenkoUpdate",PARAM_TYPE_BOOL, offsetof(estim_config_t, useBorisenkoUpdate),
   "use Borisenko et al. (2019) parameter update algorithm in algorithm EE"},

  {"learningRate",     PARAM_TYPE_DOUBLE,offsetof(estim_config_t, learningRate),
   "learning rate a in Borisenko update step of algorithm EE"},

  {"minTheta",         PARAM_TYPE_DOUBLE,offsetof(estim_config_t, minTheta),
   "min abs value of theta to stop zero in Borisenko EE algorithm update step"},

  {"computeStats",     PARAM_TYPE_BOOL,offsetof(estim_config_t, computeStats),
   "compute observed statistics corresponding to parameters being estimated"},

  {"observedStatsFilePrefix", PARAM_TYPE_STRING,
   offsetof(estim_config_t, obs_stats_file_prefix),
   "observed sufficient statistics output filename prefix"},

  {"outputFileSuffixBase", PARAM_TYPE_UINT,
   offsetof(estim_config_t, outputFileSuffixBase),
   "number to add task number to for output file suffixes"},

  {STRUCT_PARAMS_STR,  PARAM_TYPE_SET,      0, /*no offset, coded explicitly*/
  "structural parameters to estimate"},

  {ATTR_PARAMS_STR, PARAM_TYPE_SET,         0, /*no offset, coded explicitly*/
  "binary/categorical/continuous/set attribute parameters to estimate"},

  {DYADIC_PARAMS_STR, PARAM_TYPE_SET,       0, /*no offset, coded explicitly*/
  "dyadic covariate parameters to estimate"},

  {ATTR_INTERACTION_PARAMS_STR,PARAM_TYPE_SET, 0,/*no offset, coded explicitly*/
   "attribute pair interaction parameters to estimate"}
};
const uint_t NUM_ESTIM_CONFIG_PARAMS = sizeof(ESTIM_CONFIG_PARAMS) /
  sizeof(ESTIM_CONFIG_PARAMS[0]);



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
estim_config_t ESTIM_CONFIG = {
  0.1,   /* ACA_S */       
  DEFAULT_ACA_EE, /* ACA_EE */      
  DEFAULT_COMPC, /* compC */       
  1000,  /* samplerSteps */
  100,   /* Ssteps */      
  500,   /* EEsteps */     
  100,   /* EEinnerSteps */
  FALSE, /* outputAllSteps */
  FALSE, /* useIFDsampler */
  FALSE, /* useTNTsampler */
  DEFAULT_IFD_K,   /* ifd_K */
  FALSE, /* outputSimulatedNetwork */
  NULL,  /* arclist_filename */
  NULL,  /* binattr_filename */
  NULL,  /* catattr_filename */
  NULL,  /* contattr_filename */
  NULL,  /* setattr_filename */
  NULL,  /* theta_file_prefix */
  NULL,  /* dzA_file_prefix */
  NULL,  /* sim_net_file_prefix */
  NULL,  /* zone_filename */
  FALSE, /* useConditionalEstimation */
  FALSE, /* forbidReciprocity */
  FALSE, /* useBorisenkoUpdate */
  DEFAULT_LEARNING_RATE, /* learningRate */
  DEFAULT_MIN_THETA,     /* minTheta */
  FALSE, /* computeStats */
  NULL,  /* obs_stats_file_prefix */
  0,     /* outputFileSuffixBase */
  {
    0,     /* num_change_stats_funcs */
    NULL,  /* change_stats_funcs */
    NULL,  /* param_names */
    NULL,  /* param_lambdas */
    NULL,  /* param_values */
    0,     /* num_attr_change_stats_funcs */
    NULL,  /* attr_change_stats_funcs */
    NULL,  /* attr_names */
    NULL,  /* attr_indices */
    NULL,  /* attr_param_names */
    NULL,  /* attr_param_values */
    0,     /* num_dyadic_change_stats_funcs */
    NULL,  /* dyadic_change_stats_funcs */
    NULL,  /* dyadic_names */
    NULL,  /* dyadic_indices */
    NULL,  /* dyadic_types */
    NULL,  /* dyadic_param_names */
    NULL,  /* dyadic_param_values */
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
static bool ESTIM_CONFIG_IS_SET[] = {
  FALSE, /* ACA_S */       
  FALSE, /* ACA_EE */      
  FALSE, /* compC */       
  FALSE, /* samplerSteps */
  FALSE, /* Ssteps */      
  FALSE, /* EEsteps */     
  FALSE, /* EEinnerSteps */
  FALSE, /* outputAllSteps */
  FALSE, /* useIFDsampler */
  FALSE, /* useTNTsampler */
  FALSE, /* ifd_K */
  FALSE, /* outputSimulatedNetwork */
  FALSE, /* arclist_filename */
  FALSE, /* binattr_filename */
  FALSE, /* catattr_filename */
  FALSE, /* contattr_filename */
  FALSE, /* setattr_filename */
  FALSE, /* theta_file_prefix */
  FALSE, /* dzA_file_prefix */
  FALSE, /* sim_net_file_prefix */
  FALSE, /* zone_filename */
  FALSE, /* useConditionalEstimation */
  FALSE, /* forbidReciprocity */
  FALSE, /* useBorisenkoUpdate */
  FALSE, /* learningRate */
  FALSE, /* minTheta */
  FALSE, /* computeStats */
  FALSE, /* obs_stats_file_prefix */
  FALSE, /* outputFileSuffixBase */
  FALSE, /* (NOT USED) structParams */
  FALSE, /* (NOT USED) attrParams */
  FALSE, /* (NOT USED) dyadicParams */
  FALSE  /* (NOT USED) attrInteractionParams */
};
static const uint_t NUM_ESTIM_CONFIG_IS_SET = sizeof(ESTIM_CONFIG_IS_SET)/sizeof(ESTIM_CONFIG_IS_SET[0]);



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
estim_config_t *parse_estim_config_file(const char *config_filename)
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
                                    &ESTIM_CONFIG, ESTIM_CONFIG_IS_SET,
                                    &ESTIM_CONFIG.param_config,
                                    ESTIM_CONFIG_PARAMS,
                                    NUM_ESTIM_CONFIG_PARAMS, FALSE) != 0)
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
  return &ESTIM_CONFIG; /* return pointer to static CONFIG structure */
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
void free_estim_config_struct(estim_config_t *config)
{
  /* In fact parse_config_file() returns pointer to static CONFIG struct,
     so just free the pointers inside it */
  if (!config)
    return;
  assert(config == &ESTIM_CONFIG);
  free(config->arclist_filename);
  free(config->binattr_filename);
  free(config->catattr_filename);
  free(config->contattr_filename);
  free(config->setattr_filename);
  free(config->theta_file_prefix);
  free(config->dzA_file_prefix);
  free(config->sim_net_file_prefix);
  free(config->zone_filename);
  free_param_config_struct(&config->param_config);
}



/*
 * Initialize the parser. This consists only of setting the default string
 * parameter values (see comments on initialization of static CONFIG struct)
 */
void init_estim_config_parser(void)
{
  assert(NUM_ESTIM_CONFIG_IS_SET == NUM_ESTIM_CONFIG_PARAMS);
  ESTIM_CONFIG.theta_file_prefix = safe_strdup("theta_values");
  ESTIM_CONFIG.dzA_file_prefix = safe_strdup("dzA_values");
  ESTIM_CONFIG.sim_net_file_prefix = safe_strdup("sim");
  ESTIM_CONFIG.obs_stats_file_prefix = safe_strdup("obs_stats");
}



