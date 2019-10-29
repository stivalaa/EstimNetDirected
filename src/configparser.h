#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H
/*****************************************************************************
 * 
 * File:    configparser.h
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * Parse the configuration file to get algorithm parameters, input filenames,
 * parameters to estimate, etc.
 *
 * This file for constants, functions etc. shared by estimation and simulation
 * (and potentially other) configuration parsing.
 *
 * The config file is a text file with comments marked by '#'
 * character, and "keyword = value" pairs.  See config.txt for example
 * config file.
 *
 *
 ****************************************************************************/

#include <stdio.h>
#include "utils.h"
#include "changeStatisticsDirected.h"

/* These must be macros not const to use in initializer  */
/* not case sensitive */
#define STRUCT_PARAMS_STR             "structParams"
#define ATTR_PARAMS_STR               "attrParams"
#define DYADIC_PARAMS_STR             "dyadicParams"
#define ARC_PARAM_STR                 "Arc"
#define ATTR_INTERACTION_PARAMS_STR   "attrInteractionParams"



/*****************************************************************************
 *
 * constant declarations
 *
 ****************************************************************************/

extern const size_t TOKSIZE;   /* maximum size of a token */


/*****************************************************************************
 *
 * type definitions
 *
 ****************************************************************************/


/* config parameter types */
typedef enum param_type_e {
  PARAM_TYPE_INVALID,  /* invalid type, used as error return value */
  PARAM_TYPE_DOUBLE,   /* numeric (floating point) */
  PARAM_TYPE_UINT,     /* numeric (unsigned integer) */
  PARAM_TYPE_BOOL,     /* Boolean ("True" or "False" in config, bool in struct*/
  PARAM_TYPE_STRING,   /* string (may be quoted, not necessarily) */
  PARAM_TYPE_SET       /* comma delimited set of other params enclosed in {} */
} param_type_e;

/* ERGM attribute parameter type */
typedef enum attr_type_e {
  ATTR_TYPE_INVALID,        /* invalid type, used as error return value */
  ATTR_TYPE_BINARY,         /* binary attribute type (0/1)*/
  ATTR_TYPE_CATEGORICAL,    /* categorical attribute type (uint) */
  ATTR_TYPE_CONTINUOUS,     /* continuous attribute type (double) */
  ATTR_TYPE_SET             /* set attribute type (array of set_elem_e) */
} attr_type_e;

/* ERGM dyadic covariate parameter type */
typedef enum dyadic_type_e {
  DYADIC_TYPE_INVALID,       /* invalid type, used as error return value */
  DYADIC_TYPE_GEODISTANCE,   /* continuous geographic distance from lat/long */
  DYADIC_TYPE_EUCLIDEANDISTANCE /* continuous Euclidean distance from x/y/z */
} dyadic_type_e;



/* configuration parameter */
typedef struct config_param_s {
  const char   *name;  /* parameter name (keyword), not case sensitive */
  param_type_e  type;  /* parameter value type */
  size_t        offset; /* offset in config_t of field to set with value*/
  const char   *description; /* description of parameter */
} config_param_t;

/* ERGM structural parameter */
typedef struct struct_param_s {
  const char          *name;                 /* structural parameter name */
  change_stats_func_t *change_stats_func;    /* corresponding change stat */
} struct_param_t;

/* ERGM attribute parameter */
typedef struct attr_param_s {
  const char          *name;                 /* attribute parameter name */
  attr_type_e          type;                 /* attribute parameter type */
  attr_change_stats_func_t *attr_change_stats_func;  /* corresponding func. */
} attr_param_t;

/* ERGM dyadic covariate parameter */
typedef struct dyadic_param_s {
  const char          *name;                 /*dyadic covariate parameter name */
  dyadic_type_e        type;                 /*dyadic covariate parameter type */
  dyadic_change_stats_func_t *dyadic_change_stats_func; /* corresponding func. */
} dyadic_param_t;

/* ERGM attribute interaction parameter */
typedef struct attr_interaction_param_s {
  const char          *name;       /* attribute interaction parameter name */
  attr_type_e          type;       /* attribute interaction parameter type */
  attr_interaction_change_stats_func_t *attr_interaction_change_stats_func;  /* corresponding func. */
} attr_interaction_param_t;



/*
 * values built by confiparser.c functions from parsed config settings
 */
typedef struct param_config_s {
  uint_t num_change_stats_funcs;           /* length of change_stats_funcs */
  change_stats_func_t **change_stats_funcs; /* structural parameter stats */
  const char          **param_names;        /* names corresponding to above */
  double               *param_values;       /* initial values corr. to above */
  uint_t num_attr_change_stats_funcs;  /* length of attr_change_stats_funcs */
  attr_change_stats_func_t **attr_change_stats_funcs; /* attr param stats */
  char                     **attr_names; /* names of attributes for above */
  uint_t *attr_indices;   /* idx into digraph binattr/cattr/contattr for above */
  const char **attr_param_names; /* names corresponding to above two */
  uint_t num_dyadic_change_stats_funcs;  /* length of dyadic_change_stats_funcs */
  dyadic_change_stats_func_t **dyadic_change_stats_funcs;/* dyadic change stats*/
  char                       **dyadic_names; /* names corresponding to above */
  uint_t *dyadic_indices;  /* idx into digraph binattr/cattr/contattr for above */
  dyadic_type_e *dyadic_types; /* dyadic paramter type corresponding to above */
  const char **dyadic_param_names; /* names corresponding to above two */
  uint_t num_attr_interaction_change_stats_funcs;  /* length of attr_interaction_change_stats_funcs */
  attr_interaction_change_stats_func_t **attr_interaction_change_stats_funcs; /* attr interaction param stats */
  string_pair_t *attr_interaction_pair_names; /* names of pairs of attributes for above */
  uint_pair_t  *attr_interaction_pair_indices;   /* pairs of indices into digraph binattr/cattr/contattr for above */
  const char **attr_interaction_param_names; /* names corresponding to above two */
} param_config_t;




/*****************************************************************************
 *
 * function prototypes
 *
 ****************************************************************************/

int fskip(FILE *f);
int isSingleCharToken(int c);
int istokenchar(int c);
int isParamNameChar(int c);
int isValidParamname(const char *s);
char *get_token(FILE *infile, char *token);
int get_paramname_value(FILE *infile, char *paramname, char *value);

attr_type_e get_attr_param_type(const char *attrParamName);
dyadic_type_e get_dyadic_param_type(const char *dyadicParamName);
attr_type_e get_attr_interaction_param_type(const char
                                            *attrInteractionParamName);

int build_attr_indices_from_names(param_config_t *pconfig, const digraph_t *g);
int build_dyadic_indices_from_names(param_config_t *pconfig, digraph_t *g);
int build_attr_interaction_pair_indices_from_names(param_config_t *pconfig,
                                                   const digraph_t *g);

void free_param_config_struct(param_config_t *pconfig);


void dump_parameter_names(void);

void dump_config_names(const void *config,
                       const config_param_t *config_params,
                       uint_t num_config_params);

int check_and_set_param_value(const char *paramname,
                              const char *valuestr,
                              FILE *infile,
                              void *config,
                              bool *config_is_set,
                              param_config_t *pconfig,
                              const config_param_t *config_params,
                              uint_t num_config_params,
                              bool requireErgmValue);


#endif /* CONFIGPARSER_H */

