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
const config_param_t CONFIG_PARAMS[] = {
  {"ACA_S",         PARAM_TYPE_DOUBLE,   offsetof(config_t, ACA_S),
   "multiplier for step size in Algorithm S"},

  {"ACA_EE",        PARAM_TYPE_DOUBLE,   offsetof(config_t, ACA_EE),
   "multiplier for step size in Algorithm EE"},

  {"compC",         PARAM_TYPE_DOUBLE,   offsetof(config_t, compC),
   "multiplier of sd(theta)/mean(theta) to limit variance"},

  {"samplerSteps",  PARAM_TYPE_UINT,     offsetof(config_t, samplerSteps),
   "sampler iterations (per algorithm step)"},

  {"Ssteps",        PARAM_TYPE_UINT,     offsetof(config_t, Ssteps),
   "steps of Algorithm S"},

  {"EEsteps",       PARAM_TYPE_UINT,     offsetof(config_t, EEsteps),
   "steps of Algorithm EE"},

  {"EEinnerSteps",  PARAM_TYPE_UINT,     offsetof(config_t, EEinnerSteps),
   "inner iterations of Algorithm EE"},

  {"outputAllSteps", PARAM_TYPE_BOOL,    offsetof(config_t, outputAllSteps),
   "output theta and dzA values on every iteration of EE algorithm)"},

  {"useIFDsampler", PARAM_TYPE_BOOL,    offsetof(config_t, useIFDsampler),
   "use Improved Fixed Density sampler instead of basic sampler"},

  {"ifd_K",         PARAM_TYPE_DOUBLE,  offsetof(config_t, ifd_K),
   "multiplier for auxiliary parameter step size in IFD sampler"},

  {"outputSimulatedNetwork", PARAM_TYPE_BOOL,
   offsetof(config_t, outputSimulatedNetwork),
   "output simulated network in Pajek format at end of MCMC simulation"},

  {"arclistFile",   PARAM_TYPE_STRING,   offsetof(config_t, arclist_filename),
  "Network in Pajek arc list format"},

  {"binattrFile",   PARAM_TYPE_STRING,   offsetof(config_t, binattr_filename),
  "binary attributes file"},

  {"catattrFile",   PARAM_TYPE_STRING,   offsetof(config_t, catattr_filename),
  "categorical attributes file"},

  {"contattrFile",  PARAM_TYPE_STRING,   offsetof(config_t, contattr_filename),
  "continuous attributes file"},

  {"setattrFile",   PARAM_TYPE_STRING,   offsetof(config_t, setattr_filename),
  "set attributes file"},

  {"thetaFilePrefix",PARAM_TYPE_STRING,  offsetof(config_t, theta_file_prefix),
   "theta output file prefix"},

  {"dzAFilePrefix",  PARAM_TYPE_STRING,  offsetof(config_t, dzA_file_prefix),
   "dzA output file prefix"},
  
  {"simNetFilePrefix", PARAM_TYPE_STRING,
   offsetof(config_t, sim_net_file_prefix),
   "simulated network output file prefix"},

  {"zoneFile",      PARAM_TYPE_STRING,  offsetof(config_t, zone_filename),
   "snowball sample zone file"},

  {"useConditionalEstimation", PARAM_TYPE_BOOL,
   offsetof(config_t, useConditionalEstimation),
   "do conditional estimation for snowball network sample"},

  {"forbidReciprocity",PARAM_TYPE_BOOL, offsetof(config_t, forbidReciprocity),
   "constrain ERGM sampler to not allow reciprocated arcs"},
  /* This is useful for graphs that have no reciprocated arcs in the observed
     graph so cannot use Reciprocity parameter, and want to enforce 
     constraint that there can be no reciprocated arcs e.g. for a citatoin
     network. 
     TODO should have some more general way of specifying constraints
     like ergm-constraints in statnet instead of this ad-hoc way */

  {"useBorisenkoUpdate",PARAM_TYPE_BOOL, offsetof(config_t, useBorisenkoUpdate),
   "use Borisenko et al. (2019) parameter update algorithm in algorithm EE"},

  {"learningRate",     PARAM_TYPE_DOUBLE,offsetof(config_t, learningRate),
   "learning rate a in Borisenko update step of algorithm EE"},

  {"minTheta",         PARAM_TYPE_DOUBLE,offsetof(config_t, minTheta),
   "min abs value of theta to stop zero in Borisenko EE algorithm update step"},

  {STRUCT_PARAMS_STR,  PARAM_TYPE_SET,      0, /*no offset, coded explicitly*/
  "structural parameters to estimate"},

  {ATTR_PARAMS_STR, PARAM_TYPE_SET,         0, /*no offset, coded explicitly*/
  "binary/categorical/continuous/set attribute parameters to estimate"},

  {DYADIC_PARAMS_STR, PARAM_TYPE_SET,       0, /*no offset, coded explicitly*/
  "dyadic covariate parameters to estimate"},

  {ATTR_INTERACTION_PARAMS_STR,PARAM_TYPE_SET, 0,/*no offset, coded explicitly*/
   "attribute pair interaction parameters to estimate"}
};
const uint_t NUM_CONFIG_PARAMS = sizeof(CONFIG_PARAMS) /
  sizeof(CONFIG_PARAMS[0]);



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
config_t CONFIG = {
  0.1,   /* ACA_S */       
  DEFAULT_ACA_EE, /* ACA_EE */      
  DEFAULT_COMPC, /* compC */       
  1000,  /* samplerSteps */
  100,   /* Ssteps */      
  500,   /* EEsteps */     
  100,   /* EEinnerSteps */
  FALSE, /* outputAllSteps */
  FALSE, /* useIFDsampler */
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
static bool CONFIG_IS_SET[] = {
  FALSE, /* ACA_S */       
  FALSE, /* ACA_EE */      
  FALSE, /* compC */       
  FALSE, /* samplerSteps */
  FALSE, /* Ssteps */      
  FALSE, /* EEsteps */     
  FALSE, /* EEinnerSteps */
  FALSE, /* outputAllSteps */
  FALSE, /* useIFDsampler */
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
  FALSE, /* (NOT USED) structParams */
  FALSE, /* (NOT USED) attrParams */
  FALSE, /* (NOT USED) dyadicParams */
  FALSE  /* (NOT USED) attrInteractionParams */
};
static const uint_t NUM_CONFIG_IS_SET = sizeof(CONFIG_IS_SET)/sizeof(CONFIG_IS_SET[0]);



/*****************************************************************************
 *
 * local functions
 *
 ****************************************************************************/



/* 
 * Parse the structural parameters sructParams from (open read) file infile
 * These are the set type, comma delimited names of structural parameters
 * enclosed in braces, eg.:
 * "{Arc, Reciprocity, AltInStars, AltOutStars, AltKTrianglesT}"
 * The change_stats_funcs field in the CONFIG (file static) structure
 * is set to corresponding list of change statistics function pointers.
 * The STRUCT_PARAMS constant has the table of parameter names and
 * corresponding change statistic functions.
 * Return nonzero on error else zero.
 */
static int parse_struct_params(FILE *infile, param_config_t *pconfig)
{
  char        tokenbuf[TOKSIZE];
  char       *token;
  bool        found_paramname;
  bool        last_token_was_paramname = FALSE;
  uint_t      i;
  
  if (!(token = get_token(infile, tokenbuf))) {
    fprintf(stderr, "ERROR: no tokens for structParams\n");
    return 1; 
  }
  while (token && !(strlen(token) == 1 && token[0] == CLOSE_SET_CHAR)) {
    CONFIG_DEBUG_PRINT(("parse_struct_params token '%s'\n", token));
    if (last_token_was_paramname) {
      last_token_was_paramname = FALSE;
      if (strcmp(token, ",") != 0) {
        fprintf(stderr, "ERROR: structParams expecting parameter names separated by comma\n");
        return 1;
      }
    } else {
      found_paramname = FALSE;
      for (i = 0; i < NUM_STRUCT_PARAMS; i++) {
        if (strcasecmp(token, STRUCT_PARAMS[i].name) == 0) {
          found_paramname = TRUE;
          break;
        }
      }
      if (!found_paramname) {
        fprintf(stderr, "ERROR: '%s' is not a valid structural parameter "
                "name for structParams\n", token);
        return 1;
      }
      CONFIG_DEBUG_PRINT(("structParam %s\n", token));
      last_token_was_paramname = TRUE;
      pconfig->param_names = (const char **)safe_realloc(pconfig->param_names,
                                        (pconfig->num_change_stats_funcs + 1) *
                                        sizeof(const char *));
      pconfig->change_stats_funcs = (change_stats_func_t **)
        safe_realloc(pconfig->change_stats_funcs,
                  (pconfig->num_change_stats_funcs + 1) *
                     sizeof(change_stats_func_t *));
      pconfig->param_names[pconfig->num_change_stats_funcs] =
        STRUCT_PARAMS[i].name;
      pconfig->change_stats_funcs[pconfig->num_change_stats_funcs] =
        STRUCT_PARAMS[i].change_stats_func;
      pconfig->num_change_stats_funcs++;
    }
    token = get_token(infile, tokenbuf);
  }
  return 0;
}

/*
 * Parse a single attribute parameter (paramName and corresponding
 * attr_change_stats_func) in the attrParams set from infile.  These
 * are comma-delmited attribute names inside parentheses e.g.
 * "Matching(class1,class2)". See parse_attr_params() for note on note
 * validating these names yet, just setting them in
 * CONFIG.param_config.attr_names[].  Return nonzero on error else zero.
 */
static int parse_one_attr_param(const char *paramName,
                                attr_change_stats_func_t *attr_change_stats_func,
                                FILE *infile, param_config_t *pconfig)
{
  char      tokenbuf[TOKSIZE];
  char     *token;
  bool      last_token_was_attrname = FALSE;
  bool      opening = TRUE; /* true for first iteration only to expect '(' */

  if (!(token = get_token(infile, tokenbuf))) {
    fprintf(stderr, "ERROR: no tokens for attrParam %s\n", paramName);
    return 1; 
  }
  while (token && !(strlen(token) == 1 && token[0] == CLOSE_PAREN_CHAR)) {
    CONFIG_DEBUG_PRINT(("parse_one_attr_param token '%s'\n", token));
    if (opening) {
      if (!(strlen(token) == 1 && token[0] == OPEN_PAREN_CHAR)) {
        fprintf(stderr, "ERROR: expecting %c to open list of attribute "
                "names after attrParam %s but got '%s'\n",
                OPEN_PAREN_CHAR, paramName, token);
        return 1;
      }
      opening = FALSE;
    } else {
      if (last_token_was_attrname) {
        last_token_was_attrname = FALSE;
        if (strcmp(token, ",") != 0) {
          fprintf(stderr, "ERROR: attrParams %s expecting parameter names "
                  "separated by comma\n", paramName);
          return 1;
        }
      } else {
        CONFIG_DEBUG_PRINT(("attrParam %s('%s')\n", paramName, token));
        last_token_was_attrname = TRUE;
        pconfig->attr_param_names = (const char **)
          safe_realloc(pconfig->attr_param_names,
                       (pconfig->num_attr_change_stats_funcs + 1) *
                       sizeof(const char *));
        pconfig->attr_change_stats_funcs = (attr_change_stats_func_t **)
          safe_realloc(pconfig->attr_change_stats_funcs,
                       (pconfig->num_attr_change_stats_funcs + 1) *
                       sizeof(attr_change_stats_func_t *));
        pconfig->attr_names = (char **)safe_realloc(pconfig->attr_names,
              (pconfig->num_attr_change_stats_funcs + 1) * sizeof(const char *));
        pconfig->attr_param_names[pconfig->num_attr_change_stats_funcs] = paramName;
        pconfig->attr_change_stats_funcs[pconfig->num_attr_change_stats_funcs] =
          attr_change_stats_func;
        pconfig->attr_names[pconfig->num_attr_change_stats_funcs] = safe_strdup(token);
        pconfig->num_attr_change_stats_funcs++;
      }
    }
    token = get_token(infile, tokenbuf);
  }
  return 0;
  
}

/* 
 * Parse the attribute parameters attrParams from (open read) file infile.
 * These are the set type, comma delimited names of structural parameters
 * enclosed in braces, eg.:
 * "{Sender(binaryAttribute), Matching(class1,class2)}"
 *   
 * The attr_change_stats_funcs and attr_names fields in the CONFIG
 * (file static) structure is set to corresponding list of change
 * statistics function pointers.  The ATTR_PARAMS constant has the
 * table of parameter names and corresponding change statistic
 * functions.  Return nonzero on error else zero.
 *
 * Where there are multiple attributes for one effect (as in the Matching
 * example above), separate entries are created, with the same change
 * statistics function but with a different attribute name (so there is
 * always only one attribute name per change statistic function, but 
 * there can be multiple entries for the same change statstic function).
 *
 * Note that we cannot validate the attribute names yet (while parsing
 * the config file), as it requires the attributes to have been loaded
 * from the files named in the config file. So the attr_names field
 * is set, and the attr_indices are only built from the names later
 * after the attributes values have been loaded, by calling
 * build_attr_indices_from_names().
 */
static int parse_attr_params(FILE *infile, param_config_t *pconfig)
{
  char        tokenbuf[TOKSIZE];
  char       *token;
  bool        found_paramname;
  bool        last_token_was_paramname = FALSE;
  uint_t      i;
  
  if (!(token = get_token(infile, tokenbuf))) {
    fprintf(stderr, "ERROR: no tokens for attrParams\n");
    return 1; 
  }
  while (token && !(strlen(token) == 1 && token[0] == CLOSE_SET_CHAR)) {
    CONFIG_DEBUG_PRINT(("parse_attr_params token '%s'\n", token));
    if (last_token_was_paramname) {
      last_token_was_paramname = FALSE;
      if (strcmp(token, ",") != 0) {
        fprintf(stderr, "ERROR: attrParams expecting parameter names separated by comma\n");
        return 1;
      }
    } else {
      found_paramname = FALSE;
      for (i = 0; i < NUM_ATTR_PARAMS; i++) {
        if (strcasecmp(token, ATTR_PARAMS[i].name) == 0) {
          found_paramname = TRUE;
          break;
        }
      }
      if (!found_paramname) {
        fprintf(stderr, "ERROR: '%s' is not a valid attribute parameter "
                "name for attrParams\n", token);
        return 1;
      }
      CONFIG_DEBUG_PRINT(("attrParam %s\n", token));
      last_token_was_paramname = TRUE;
      if (parse_one_attr_param(ATTR_PARAMS[i].name,
                               ATTR_PARAMS[i].attr_change_stats_func,
                               infile, pconfig)) {
        fprintf(stderr, "ERROR parsing attrParams %s\n", ATTR_PARAMS[i].name);
        return 1;
      }
          }
    token = get_token(infile, tokenbuf);
  }
  return 0;
}




/*
 * Parse a single dyadic covariate parameter (paramName and corresponding
 * attr_change_stats_func) in the dyadicParams set from infile.  These
 * are comma-delmited attribute names inside parentheses e.g.
 * Note that we cannot validate the attribute names yet (while parsing
 * the config file), as it requires the attributes to have been loaded
 * from the files named in the config file. So the attr_names field
 * is set, and the attr_indices are only built from the names later
 * after the attributes values have been loaded, by calling
 * build_dyadic_indices_from_names().
 * Return nonzero on error else zero.
 */
static int parse_one_dyadic_param(const char *paramName,
                                dyadic_change_stats_func_t *dyadic_change_stats_func,
                                  FILE *infile, param_config_t *pconfig)
{
  char      tokenbuf[TOKSIZE];
  char     *token;
  bool      last_token_was_attrname = FALSE;
  bool      opening = TRUE; /* true for first iteration only to expect '(' */

  if (!(token = get_token(infile, tokenbuf))) {
    fprintf(stderr, "ERROR: no tokens for dyadicParam %s\n", paramName);
    return 1; 
  }
  while (token && !(strlen(token) == 1 && token[0] == CLOSE_PAREN_CHAR)) {
    CONFIG_DEBUG_PRINT(("parse_one_dyadic_param token '%s'\n", token));
    if (opening) {
      if (!(strlen(token) == 1 && token[0] == OPEN_PAREN_CHAR)) {
        fprintf(stderr, "ERROR: expecting %c to open list of attribute "
                "names after dyadicParam %s but got '%s'\n",
                OPEN_PAREN_CHAR, paramName, token);
        return 1;
      }
      opening = FALSE;
    } else {
      if (last_token_was_attrname) {
        last_token_was_attrname = FALSE;
        if (strcmp(token, ",") != 0) {
          fprintf(stderr, "ERROR: dyadicParams %s expecting parameter names "
                  "separated by comma\n", paramName);
          return 1;
        }
      } else {
        CONFIG_DEBUG_PRINT(("dyadicParam %s('%s')\n", paramName, token));
        last_token_was_attrname = TRUE;
        pconfig->dyadic_param_names = (const char **)
          safe_realloc(pconfig->dyadic_param_names,
                       (pconfig->num_dyadic_change_stats_funcs + 1) *
                       sizeof(const char *));
        pconfig->dyadic_change_stats_funcs = (dyadic_change_stats_func_t **)
          safe_realloc(pconfig->dyadic_change_stats_funcs,
                       (pconfig->num_dyadic_change_stats_funcs + 1) *
                       sizeof(dyadic_change_stats_func_t *));
        pconfig->dyadic_names = (char **)safe_realloc(pconfig->dyadic_names,
              (pconfig->num_dyadic_change_stats_funcs + 1) * sizeof(char *));
        pconfig->dyadic_param_names[pconfig->num_dyadic_change_stats_funcs] = paramName;
        pconfig->dyadic_change_stats_funcs[pconfig->num_dyadic_change_stats_funcs] =
          dyadic_change_stats_func;
        pconfig->dyadic_names[pconfig->num_dyadic_change_stats_funcs] = safe_strdup(token);
        pconfig->num_dyadic_change_stats_funcs++;
      }
    }
    token = get_token(infile, tokenbuf);
  }
  return 0;
  
}

/* 
 * Parse the attribute parameters dyadicParams from (open read) file infile.
 * These are the set type, comma delimited names of structural parameters
 * enclosed in braces, eg.:
 *   
 * The dyadic_change_stats_funcs and dyadic_names fields in the CONFIG
 * (file static) structure is set to corresponding list of change
 * statistics function pointers.  The DYADIC_PARAMS constant has the
 * table of parameter names and corresponding change statistic
 * functions.  Return nonzero on error else zero.
 *
 */
static int parse_dyadic_params(FILE *infile, param_config_t *pconfig)
{
  char        tokenbuf[TOKSIZE];
  char       *token;
  bool        found_paramname;
  bool        last_token_was_paramname = FALSE;
  uint_t      i;
  
  if (!(token = get_token(infile, tokenbuf))) {
    fprintf(stderr, "ERROR: no tokens for dyadicParams\n");
    return 1; 
  }
  while (token && !(strlen(token) == 1 && token[0] == CLOSE_SET_CHAR)) {
    CONFIG_DEBUG_PRINT(("parse_dyadic_params token '%s'\n", token));
    if (last_token_was_paramname) {
      last_token_was_paramname = FALSE;
      if (strcmp(token, ",") != 0) {
        fprintf(stderr, "ERROR: dyadicParams expecting parameter names separated by comma\n");
        return 1;
      }
    } else {
      found_paramname = FALSE;
      for (i = 0; i < NUM_DYADIC_PARAMS; i++) {
        if (strcasecmp(token, DYADIC_PARAMS[i].name) == 0) {
          found_paramname = TRUE;
          break;
        }
      }
      if (!found_paramname) {
        fprintf(stderr, "ERROR: '%s' is not a valid attribute parameter "
                "name for dyadicParams\n", token);
        return 1;
      }
      CONFIG_DEBUG_PRINT(("dyadicParam %s\n", token));
      last_token_was_paramname = TRUE;
      if (parse_one_dyadic_param(DYADIC_PARAMS[i].name,
                                 DYADIC_PARAMS[i].dyadic_change_stats_func,
                                 infile, pconfig)) {
        fprintf(stderr, "ERROR parsing dyadicParams %s\n", DYADIC_PARAMS[i].name);
        return 1;
      }
          }
    token = get_token(infile, tokenbuf);
  }
  return 0;
}




/*
 * Parse a single attrInteractionParams attribute interaction
 * parameter (paramName and corresponding attr_change_stats_func) in
 * the attrParams set from infile.  These are a single comma-delmited
 * paior of attribute names inside parentheses e.g.
 * "MatchingInteraction(class1,class2)".  See
 * parse_attr_interaction_params() for note on note validating these
 * names yet, just setting them in
 * CONFIG.param_config.attr_interaction_pair_names[].  Return nonzero on error else
 * zero.
 */
static int parse_one_attr_interaction_param(const char *paramName,
                  attr_interaction_change_stats_func_t *attr_interaction_change_stats_func,
                  FILE *infile, param_config_t *pconfig)
{
  char      tokenbuf[TOKSIZE];
  char     *token;
  bool      last_token_was_attrname = FALSE;
  bool      opening = TRUE; /* true for first iteration only to expect '(' */
  uint_t    num_attr_names = 0;

  if (!(token = get_token(infile, tokenbuf))) {
    fprintf(stderr, "ERROR: no tokens for attrInteractionParam %s\n",
            paramName);
    return 1; 
  }
  while (token && !(strlen(token) == 1 && token[0] == CLOSE_PAREN_CHAR)) {
    CONFIG_DEBUG_PRINT(("parse_one_attr_interaction_param token '%s'\n", token));
    if (opening) {
      if (!(strlen(token) == 1 && token[0] == OPEN_PAREN_CHAR)) {
        fprintf(stderr, "ERROR: expecting %c to open list of attribute "
                "names after attrParam %s but got '%s'\n",
                OPEN_PAREN_CHAR, paramName, token);
        return 1;
      }
      opening = FALSE;
    } else {
      if (last_token_was_attrname) {
        last_token_was_attrname = FALSE;
        if (strcmp(token, ",") != 0) {
          fprintf(stderr, "ERROR: attrInteractionParams %s expecting "
                  "two parameter names "
                  "separated by comma\n", paramName);
          return 1;
        }
      } else {
        CONFIG_DEBUG_PRINT(("attrInteractionParam %s('%s') [%u]\n",
                            paramName, token, num_attr_names));
        last_token_was_attrname = TRUE;
        if (num_attr_names == 0) {
          pconfig->attr_interaction_param_names = (const char **)
            safe_realloc(pconfig->attr_interaction_param_names,
                         (pconfig->num_attr_interaction_change_stats_funcs + 1) *
                         sizeof(const char *));
          pconfig->attr_interaction_change_stats_funcs =
            (attr_interaction_change_stats_func_t **)
            safe_realloc(pconfig->attr_interaction_change_stats_funcs,
                         (pconfig->num_attr_interaction_change_stats_funcs + 1) *
                         sizeof(attr_interaction_change_stats_func_t *));
          pconfig->attr_interaction_pair_names = (string_pair_t *)safe_realloc(
            pconfig->attr_interaction_pair_names,
            (pconfig->num_attr_interaction_change_stats_funcs + 1) * sizeof(string_pair_t));
          pconfig->attr_interaction_pair_names[
            pconfig->num_attr_interaction_change_stats_funcs].first = safe_strdup(token);
          pconfig->attr_interaction_pair_names[
            pconfig->num_attr_interaction_change_stats_funcs].second = NULL;
          pconfig->attr_interaction_param_names[pconfig->num_attr_interaction_change_stats_funcs] = paramName;
          pconfig->attr_interaction_change_stats_funcs[pconfig->num_attr_interaction_change_stats_funcs] =
            attr_interaction_change_stats_func;
          num_attr_names++;
        }  else if (num_attr_names == 1) {
          pconfig->attr_interaction_pair_names[
            pconfig->num_attr_interaction_change_stats_funcs].second = safe_strdup(token);
          num_attr_names++;
          pconfig->num_attr_interaction_change_stats_funcs++;
        } else {
          fprintf(stderr, "ERROR: attrInteractionParams %s expecting "
                  "exactly two parameter names "
                  "separated by comma\n", paramName);
          return 1;
        }
      }
    }
    token = get_token(infile, tokenbuf);
  }
  if (num_attr_names != 2) {
    fprintf(stderr, "ERROR: attrInteractionParams %s was expecting "
            "exactly two parameter names separated by a comma\n", paramName);
    return 1;
  }
  return 0;
}


/* 
 * Parse the attribute parameters attrInteractionParams from (open
 * read) file infile. These are the set type, comma delimited names
 * of parameters enclosed in braces, eg.:
 * "{MatchingInteraction(class1,class2)}"
 *   
 * The attr_interaction_change_stats_funcs and
 * attr_intearction_pair_names fields in the CONFIG (file static)
 * structure is set to corresponding list of change statistics
 * function pointers.  The ATTR_INTERACTION_PARAMS constant has the
 * table of parameter names and corresponding change statistic
 * functions.  Return nonzero on error else zero.
 *
 * Note that we cannot validate the attribute names yet (while parsing
 * the config file), as it requires the attributes to have been loaded
 * from the files named in the config file. So the
 * attr_interaction_pairnames field is set, and the attr_indices are
 * only built from the names later after the attributes values have
 * been loaded, by calling
 * build_attr_interaction_indices_from_names().
 */
static int parse_attr_interaction_params(FILE *infile, param_config_t *pconfig)
{
  char        tokenbuf[TOKSIZE];
  char       *token;
  bool        found_paramname;
  bool        last_token_was_paramname = FALSE;
  uint_t      i;
  
  if (!(token = get_token(infile, tokenbuf))) {
    fprintf(stderr, "ERROR: no tokens for attrInteractionParams\n");
    return 1; 
  }
  while (token && !(strlen(token) == 1 && token[0] == CLOSE_SET_CHAR)) {
    CONFIG_DEBUG_PRINT(("parse_attr_interaction_interaction_params "
                        "token '%s'\n", token));
    if (last_token_was_paramname) {
      last_token_was_paramname = FALSE;
      if (strcmp(token, ",") != 0) {
        fprintf(stderr, "ERROR: attrIntearctionParams expecting "
                "parameter names separated by comma\n");
        return 1;
      }
    } else {
      found_paramname = FALSE;
      for (i = 0; i < NUM_ATTR_INTERACTION_PARAMS; i++) {
        if (strcasecmp(token, ATTR_INTERACTION_PARAMS[i].name) == 0) {
          found_paramname = TRUE;
          break;
        }
      }
      if (!found_paramname) {
        fprintf(stderr, "ERROR: '%s' is not a valid attribute parameter "
                "name for attrInteractionParams\n", token);
        return 1;
      }
      CONFIG_DEBUG_PRINT(("attrInteractionParam %s\n", token));
      last_token_was_paramname = TRUE;
      if (parse_one_attr_interaction_param(ATTR_INTERACTION_PARAMS[i].name,
                               ATTR_INTERACTION_PARAMS[i].attr_interaction_change_stats_func,
                                           infile, pconfig)) {
        fprintf(stderr, "ERROR parsing attrInteractionParams %s\n",
                ATTR_INTERACTION_PARAMS[i].name);
        return 1;
      }
          }
    token = get_token(infile, tokenbuf);
  }
  return 0;
}




/* 
 * Given the parameter name and value parsed from the config file,
 * check that they are valid, using the file global constant
 * CONFIG_PARAMS and set the corresponding fields in the file static
 * CONFIG structure. Returns nonzero on error, else zero.
 * The infile (open read) from which they are parsed is also
 * passed.  For set values, the valuestr is just '{' and more parsing
 * is required, for which the infile is passed to the appropriate function.
 */
 
static int check_and_set_param_value(const char *paramname,
                                     const char *valuestr,
                                     FILE *infile)
{
  uint_t i;
  bool   found_paramname = FALSE;
  double valuedouble;
  char  *endptr; /* for strtod() */
  uint_t valueint;
  bool   valuebool;

  for (i = 0; i < NUM_CONFIG_PARAMS; i++) {
    if (strcasecmp(paramname, CONFIG_PARAMS[i].name) == 0) {
      found_paramname = TRUE;
      break;
    }
  }
  if (!found_paramname) {
    fprintf(stderr, "ERROR: invalid parameter name '%s'\n", paramname);
    return 1;
  }

  if (CONFIG_IS_SET[i]) {
    fprintf(stderr, "ERROR: parameter %s is set more than once\n",
            CONFIG_PARAMS[i].name);
    return 1;
  }
  CONFIG_IS_SET[i] = TRUE;
  
  switch(CONFIG_PARAMS[i].type) {
    case PARAM_TYPE_DOUBLE:  /* numeric (floating point) */
      valuedouble = strtod(valuestr, &endptr);
      if (*endptr != '\0') {
        fprintf(stderr, "ERROR: expecting floating point value for parameter %s but got '%s'\n", paramname, valuestr);
        return 1;
      }
      CONFIG_DEBUG_PRINT(("%s = %lg\n", paramname, valuedouble));
      *(double *)((char *)&CONFIG + CONFIG_PARAMS[i].offset) = valuedouble;
      break;
      
    case PARAM_TYPE_UINT:    /* numeric (unsigned integer) */
      if (sscanf(valuestr, "%u", &valueint) != 1) {
        fprintf(stderr, "ERROR: expecting unsigned integer value for parameter %s but got '%s'\n", paramname, valuestr);
        return 1;
      }
      CONFIG_DEBUG_PRINT(("%s = %u\n", paramname, valueint));
      *(uint_t *)((char *)&CONFIG + CONFIG_PARAMS[i].offset) = valueint;
      break;

    case PARAM_TYPE_BOOL:    /* Boolean */
      if (strcasecmp(valuestr, TRUE_STR) == 0)
        valuebool = TRUE;
      else if (strcasecmp(valuestr, FALSE_STR) == 0)
        valuebool = FALSE;
      else {
        fprintf(stderr, "ERROR: expecting Boolean value for parameter %s but got '%s'\n", paramname, valuestr);
        return 1;
      }
      CONFIG_DEBUG_PRINT(("%s = %d\n", paramname, valuebool));
      *(bool *)((char *)&CONFIG + CONFIG_PARAMS[i].offset) = valuebool;
      break;
      
    case PARAM_TYPE_STRING:  /* string (may be quoted, not necessarily) */
      /* free default/previous value first (may be NULL, which is alright)*/
      free(*(char **)((char *)&CONFIG + CONFIG_PARAMS[i].offset));
      *(char **)((char *)&CONFIG + CONFIG_PARAMS[i].offset) =
        safe_strdup(valuestr);
      CONFIG_DEBUG_PRINT(("%s = %s\n", paramname, valuestr));
      /* TODO handle quoted string */
      break;
      
    case PARAM_TYPE_SET:
      /* comma delimited set of other params enclosed in {}, so rather
         than using valuestr here, which is just empty, parse the tokens
         in the set
      */
      if (strcasecmp(paramname, STRUCT_PARAMS_STR) == 0) {
        if (strlen(valuestr) != 1 || valuestr[0] != OPEN_SET_CHAR) {
          fprintf(stderr,
                  "ERROR: expecting %c for structParams but got '%s'\n",
                  OPEN_SET_CHAR, valuestr);
          return 1;
        }
        if (CONFIG.param_config.num_change_stats_funcs > 0) {
          fprintf(stderr, "ERROR: %s specified more than once\n",
                  STRUCT_PARAMS_STR);
          return 1;
        } 
        return parse_struct_params(infile, &CONFIG.param_config); 
      } else if (strcasecmp(paramname, ATTR_PARAMS_STR) == 0) {
        if (CONFIG.param_config.num_attr_change_stats_funcs > 0) {
          fprintf(stderr, "ERROR: %s specified more than once\n",
                  ATTR_PARAMS_STR);
          return 1;
        }
        return parse_attr_params(infile, &CONFIG.param_config);
      } else if (strcasecmp(paramname, DYADIC_PARAMS_STR) == 0) {        
        if (CONFIG.param_config.num_dyadic_change_stats_funcs > 0) {
          fprintf(stderr, "ERROR: %s specified more than once\n",
                  DYADIC_PARAMS_STR);
          return 1;
        }
        return parse_dyadic_params(infile, &CONFIG.param_config);
      } else if (strcasecmp(paramname, ATTR_INTERACTION_PARAMS_STR) == 0) {
        if (CONFIG.param_config.num_attr_interaction_change_stats_funcs > 0) {
          fprintf(stderr, "ERROR: %s specified more than once\n",
                  ATTR_INTERACTION_PARAMS_STR);
          return 1;
        }
        return parse_attr_interaction_params(infile, &CONFIG.param_config);
      } else {
        fprintf(stderr, "ERROR (internal): unknown parameter %s\n", paramname);
        return 1;
      }
      break;
        
    default:
      fprintf(stderr, "ERROR (internal): unknown parameter type %d\n",
              CONFIG_PARAMS[i].type);
      return 1;
      break;
  }
  return 0;
}


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
config_t *parse_config_file(const char *config_filename)
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
      if (check_and_set_param_value(paramname, value, config_file) != 0)
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
  return &CONFIG; /* return pointer to static CONFIG structure */
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
void free_config_struct(config_t *config)
{
  /* In fact parse_config_file() returns pointer to static CONFIG struct,
     so just free the pointers inside it */
  if (!config)
    return;
  assert(config == &CONFIG);
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
void init_config_parser(void)
{
  assert(NUM_CONFIG_IS_SET == NUM_CONFIG_PARAMS);
  CONFIG.theta_file_prefix = safe_strdup("theta_values");
  CONFIG.dzA_file_prefix = safe_strdup("dzA_values");
  CONFIG.sim_net_file_prefix = safe_strdup("sim");
}



