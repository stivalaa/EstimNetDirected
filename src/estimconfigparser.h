#ifndef ESTIMCONFIGPARSER_H
#define ESTIMCONFIGPARSER_H
/*****************************************************************************
 * 
 * File:    estimconfigparser.h
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * Parse the configuration file to get algorithm parameters, input filenames,
 * parameters to estimate, etc.
 *
 * The config file is a text file with comments marked by '#'
 * character, and "keyword = value" pairs.  See config.txt for example
 * config file.
 *
 *
 ****************************************************************************/

#include "configparser.h"


/* These must be macros not const to use in initializer  */

#define DEFAULT_ACA_EE        1e-09   /* default value for ACA_EE */
#define DEFAULT_COMPC         1e-02   /* default value for compC */
#define DEFAULT_LEARNING_RATE 0.001   /* default value of learningRate */
#define DEFAULT_MIN_THETA     0.01    /* default value of minTheta */


/*****************************************************************************
 *
 * type definitions
 *
 ****************************************************************************/

typedef struct estim_config_s {
  /*
   * Parameters parsed directly from config file
   */
  double ACA_S;           /* multiplier for step size in Algorithm S */
  double ACA_EE;          /* multiplier for step size in Algorithm EE */
  double compC;           /* multiplier of sd/mean theta to limit variance */
  uint_t samplerSteps;    /* sampler iterations per algorithm step */
  uint_t Ssteps;          /* steps of Algorithm S */
  uint_t EEsteps;         /* steps of Algorithm EE */
  uint_t EEinnerSteps;    /* inner iterations of Algorithm EE */
  bool   outputAllSteps;   /* write theta and dzA every iteration not just outer*/
  bool   useIFDsampler;   /* Use IFD sampler instead of basic sampler */
  bool   useTNTsampler;   /* Use TNT sampler (not basic or IFD sampler) */
  double ifd_K;           /* multiplier for aux parameter step size in IFD sampler */
  bool  outputSimulatedNetwork; /* output simulated network at end */
  char *arclist_filename; /* filename of Pajek file with digraph to estimate */
  char *binattr_filename; /* filename of binary attributes file or NULL */
  char *catattr_filename; /* filename of categorical attributes file or NULL */
  char *contattr_filename;/* filename of continuous attributes file or NULL */
  char *setattr_filename; /* filename of set attributes file or NULL */
  char *theta_file_prefix;/* theta output filename prefix */
  char *dzA_file_prefix;  /* dzA output filename prefix */
  char *sim_net_file_prefix; /* simulated network output filename prefix */
  char *zone_filename;    /* filename of snowball sampling zone file or NULL */
  bool  useConditionalEstimation; /*conditional estimation of snowball sample */ 
  bool  forbidReciprocity; /* do not allow reciprocated arcs in sampler */
  bool  useBorisenkoUpdate; /* use Borisenko et al. update algorithm */
  double learningRate;      /* learning rate (multiplier) in Borisenko update */
  double minTheta;          /* minimum abs theta value in Borisenko update */
  bool  computeStats;       /* compute observed statistics in digraph */
  char *obs_stats_file_prefix; /* observed stats output filename prefix */
  uint_t outputFileSuffixBase;/* task number added to this for output suffixes*/
  char *term_filename;     /* filename of citation ERGM term file or NULL */
  bool  citationERGM;      /* use cERGM conditional estimation on terms */
  /*
   * values built by confiparser.c functions from parsed config settings
   */
  param_config_t param_config;
} estim_config_t;


/*****************************************************************************
 *
 * constant declarations
 *
 ****************************************************************************/

extern const config_param_t ESTIM_CONFIG_PARAMS[];
extern const uint_t NUM_ESTIM_CONFIG_PARAMS;

/*****************************************************************************
 *
 * externally visible variable declarations
 *
 ****************************************************************************/

extern estim_config_t ESTIM_CONFIG;

/*****************************************************************************
 *
 * function prototypes
 *
 ****************************************************************************/

estim_config_t *parse_estim_config_file(const char *config_filename);

void free_estim_config_struct(estim_config_t *config);

void init_estim_config_parser(void);



#endif /* ESTIMCONFIGPARSER_H */
