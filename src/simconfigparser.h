#ifndef SIMCONFIGPARSER_H
#define SIMCONFIGPARSER_H
/*****************************************************************************
 * 
 * File:    simconfigparser.h
 * Author:  Alex Stivala
 * Created: October 2019
 *
 * Parse the simulation configuration file to get algorithm
 * parameters, input filenames, parameters to estimate, etc.
 *
 * The config file is a text file with comments marked by '#'
 * character, and "keyword = value" pairs.  See config.txt for example
 * config file.
 *
 *
 ****************************************************************************/

#include "configparser.h"


/* These must be macros not const to use in initializer  */

#define SIM_DEFAULT_IFD_K         0.1     /* default value of ifd_K  */
#define SIM_DEFAULT_SAMPLE_SIZE   1000    /* sampleSize */
#define SIM_DEFAULT_INTERVAL      1000    /* interval */
#define SIM_DEFAULT_BURNIN        1000    /* burnin */

/*****************************************************************************
 *
 * type definitions
 *
 ****************************************************************************/

typedef struct sim_config_s {
  /*
   * Parameters parsed directly from config file
   */

  uint_t numNodes;        /* number of nodes in digraph */
  uint_t sampleSize;      /* number of network samples */
  uint_t interval;        /* interval (iterations) between samples */
  uint_t burnin;          /* iterations to throw out before 1st sample */
  bool   useIFDsampler;   /* Use IFD sampler instead of basic sampler */
  bool   useTNTsampler;   /* Use TNT sampler (not basic or IFD sampler) */
  double ifd_K;           /* multiplier for aux parameter step size in IFD sampler */
  bool  outputSimulatedNetworks; /* output simulated networks  */
  char *binattr_filename; /* filename of binary attributes file or NULL */
  char *catattr_filename; /* filename of categorical attributes file or NULL */
  char *contattr_filename;/* filename of continuous attributes file or NULL */
  char *setattr_filename; /* filename of set attributes file or NULL */
  char *stats_filename;   /* statistics output filename */
  char *sim_net_file_prefix; /* simulated network output filename prefix */
  char *zone_filename;    /* filename of snowball sampling zone file or NULL */
  bool  useConditionalSimulation; /*conditional simulation of snowball sample */
  bool  forbidReciprocity; /* do not allow reciprocated arcs in sampler */
  uint_t numArcs;         /* number of arcs for IFD simulation (fixed density)*/
  char *term_filename;     /* filename of citation ERGM term file or NULL */
  bool  citationERGM;      /* use cERGM conditional estimation on terms */
  char *arclist_filename;  /* filename of Pajek file for citationERGM */  
  bool allowLoops;         /* allow self-edges (loops) */
  bool isDirected;         /* directed graph (else undirected) */  
  
  /*
   * values built by confiparser.c functions from parsed config settings
   */
  param_config_t param_config;
} sim_config_t;


/*****************************************************************************
 *
 * constant declarations
 *
 ****************************************************************************/

extern const config_param_t SIM_CONFIG_PARAMS[];
extern const uint_t NUM_SIM_CONFIG_PARAMS;

/*****************************************************************************
 *
 * externally visible variable declarations
 *
 ****************************************************************************/

extern sim_config_t SIM_CONFIG;

/*****************************************************************************
 *
 * function prototypes
 *
 ****************************************************************************/

sim_config_t *parse_sim_config_file(const char *config_filename);

void free_sim_config_struct(sim_config_t *config);

void init_sim_config_parser(void);



#endif /* SIMCONFIGPARSER_H */
