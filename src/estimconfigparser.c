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
 * constants
 *
 ****************************************************************************/


/* 
 * Configuration parameter names, types and descriptions. These 
 * define what is parsed from the config file, and the help/error messages.
 * Configuration parameter names (keywords) are not case sensitive.
 */
static const config_param_t CONFIG_PARAMS[] = {
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
static const uint_t NUM_CONFIG_PARAMS = sizeof(CONFIG_PARAMS) /
  sizeof(CONFIG_PARAMS[0]);



/*****************************************************************************
 *
 * file static variables
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
};


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
bool CONFIG_IS_SET[] = {
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
static int parse_struct_params(FILE *infile)
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
      CONFIG.param_names = (const char **)safe_realloc(CONFIG.param_names,
                                        (CONFIG.num_change_stats_funcs + 1) *
                                        sizeof(const char *));
      CONFIG.change_stats_funcs = (change_stats_func_t **)
        safe_realloc(CONFIG.change_stats_funcs,
                  (CONFIG.num_change_stats_funcs + 1) *
                     sizeof(change_stats_func_t *));
      CONFIG.param_names[CONFIG.num_change_stats_funcs] =
        STRUCT_PARAMS[i].name;
      CONFIG.change_stats_funcs[CONFIG.num_change_stats_funcs] =
        STRUCT_PARAMS[i].change_stats_func;
      CONFIG.num_change_stats_funcs++;
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
 * CONFIG.attr_names[].  Return nonzero on error else zero.
 */
static int parse_one_attr_param(const char *paramName,
                                attr_change_stats_func_t *attr_change_stats_func,
                                FILE *infile)
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
        CONFIG.attr_param_names = (const char **)
          safe_realloc(CONFIG.attr_param_names,
                       (CONFIG.num_attr_change_stats_funcs + 1) *
                       sizeof(const char *));
        CONFIG.attr_change_stats_funcs = (attr_change_stats_func_t **)
          safe_realloc(CONFIG.attr_change_stats_funcs,
                       (CONFIG.num_attr_change_stats_funcs + 1) *
                       sizeof(attr_change_stats_func_t *));
        CONFIG.attr_names = (char **)safe_realloc(CONFIG.attr_names,
              (CONFIG.num_attr_change_stats_funcs + 1) * sizeof(const char *));
        CONFIG.attr_param_names[CONFIG.num_attr_change_stats_funcs] = paramName;
        CONFIG.attr_change_stats_funcs[CONFIG.num_attr_change_stats_funcs] =
          attr_change_stats_func;
        CONFIG.attr_names[CONFIG.num_attr_change_stats_funcs] = safe_strdup(token);
        CONFIG.num_attr_change_stats_funcs++;
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
static int parse_attr_params(FILE *infile)
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
                               infile)) {
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
                                FILE *infile)
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
        CONFIG.dyadic_param_names = (const char **)
          safe_realloc(CONFIG.dyadic_param_names,
                       (CONFIG.num_dyadic_change_stats_funcs + 1) *
                       sizeof(const char *));
        CONFIG.dyadic_change_stats_funcs = (dyadic_change_stats_func_t **)
          safe_realloc(CONFIG.dyadic_change_stats_funcs,
                       (CONFIG.num_dyadic_change_stats_funcs + 1) *
                       sizeof(dyadic_change_stats_func_t *));
        CONFIG.dyadic_names = (char **)safe_realloc(CONFIG.dyadic_names,
              (CONFIG.num_dyadic_change_stats_funcs + 1) * sizeof(char *));
        CONFIG.dyadic_param_names[CONFIG.num_dyadic_change_stats_funcs] = paramName;
        CONFIG.dyadic_change_stats_funcs[CONFIG.num_dyadic_change_stats_funcs] =
          dyadic_change_stats_func;
        CONFIG.dyadic_names[CONFIG.num_dyadic_change_stats_funcs] = safe_strdup(token);
        CONFIG.num_dyadic_change_stats_funcs++;
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
static int parse_dyadic_params(FILE *infile)
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
                               infile)) {
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
 * CONFIG.attr_interaction_pair_names[].  Return nonzero on error else
 * zero.
 */
static int parse_one_attr_interaction_param(const char *paramName,
                  attr_interaction_change_stats_func_t *attr_interaction_change_stats_func,
                  FILE *infile)
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
          CONFIG.attr_interaction_param_names = (const char **)
            safe_realloc(CONFIG.attr_interaction_param_names,
                         (CONFIG.num_attr_interaction_change_stats_funcs + 1) *
                         sizeof(const char *));
          CONFIG.attr_interaction_change_stats_funcs =
            (attr_interaction_change_stats_func_t **)
            safe_realloc(CONFIG.attr_interaction_change_stats_funcs,
                         (CONFIG.num_attr_interaction_change_stats_funcs + 1) *
                         sizeof(attr_interaction_change_stats_func_t *));
          CONFIG.attr_interaction_pair_names = (string_pair_t *)safe_realloc(
            CONFIG.attr_interaction_pair_names,
            (CONFIG.num_attr_interaction_change_stats_funcs + 1) * sizeof(string_pair_t));
          CONFIG.attr_interaction_pair_names[
            CONFIG.num_attr_interaction_change_stats_funcs].first = safe_strdup(token);
          CONFIG.attr_interaction_pair_names[
            CONFIG.num_attr_interaction_change_stats_funcs].second = NULL;
          CONFIG.attr_interaction_param_names[CONFIG.num_attr_interaction_change_stats_funcs] = paramName;
          CONFIG.attr_interaction_change_stats_funcs[CONFIG.num_attr_interaction_change_stats_funcs] =
            attr_interaction_change_stats_func;
          num_attr_names++;
        }  else if (num_attr_names == 1) {
          CONFIG.attr_interaction_pair_names[
            CONFIG.num_attr_interaction_change_stats_funcs].second = safe_strdup(token);
          num_attr_names++;
          CONFIG.num_attr_interaction_change_stats_funcs++;
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
static int parse_attr_interaction_params(FILE *infile)
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
                               infile)) {
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
        if (CONFIG.num_change_stats_funcs > 0) {
          fprintf(stderr, "ERROR: %s specified more than once\n",
                  STRUCT_PARAMS_STR);
          return 1;
        } 
        return parse_struct_params(infile); 
      } else if (strcasecmp(paramname, ATTR_PARAMS_STR) == 0) {
        if (CONFIG.num_attr_change_stats_funcs > 0) {
          fprintf(stderr, "ERROR: %s specified more than once\n",
                  ATTR_PARAMS_STR);
          return 1;
        }
        return parse_attr_params(infile);
      } else if (strcasecmp(paramname, DYADIC_PARAMS_STR) == 0) {        
        if (CONFIG.num_dyadic_change_stats_funcs > 0) {
          fprintf(stderr, "ERROR: %s specified more than once\n",
                  DYADIC_PARAMS_STR);
          return 1;
        }
        return parse_dyadic_params(infile);
      } else if (strcasecmp(paramname, ATTR_INTERACTION_PARAMS_STR) == 0) {
        if (CONFIG.num_attr_interaction_change_stats_funcs > 0) {
          fprintf(stderr, "ERROR: %s specified more than once\n",
                  ATTR_INTERACTION_PARAMS_STR);
          return 1;
        }
        return parse_attr_interaction_params(infile);
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
 * build_attr_indices_from_names() is called after the config file
 * is parsed by parse_config_file() and also after the 
 * attributes are loaded by calling load_digraph_from_arclist_file()
 * since the relevant filenames have to be parsed from the config
 * file before their contents can be loaded. (Could enforce that
 * those filenames are earlier in the config file than the attribute
 * parameters and do it all in one pass, but did not do so).
 *
 * The attribute names specified for attribute parameter estimation
 * are searched for in the names loaded from the attribute data files
 * stored in the digraph structure g, and if they are present, the
 * indices stored in the config.attr_indices array. Returns nonzero on
 * error else 0.
 *
 * Note that the types of attributes (binary, categorical,continuous,
 * set) are mixed together in the config.attr_indices[] array; the
 * index could be into g->binattr or g->catattr or g->contattr or
 * g->setattr, but because each change statistics function is for one
 * particular type, it uses the index in the correct array without
 * ambiguity.
 */
int build_attr_indices_from_names(config_t *config, const digraph_t *g)
{
  uint_t i, j;
  bool   found;

  config->attr_indices = safe_malloc(config->num_attr_change_stats_funcs *
                                     sizeof(uint_t));
  
  for (i = 0; i < config->num_attr_change_stats_funcs; i++) {
    found = FALSE;
    switch (get_attr_param_type(config->attr_param_names[i])) {
      case ATTR_TYPE_BINARY:
        for (j = 0; j < g->num_binattr; j++) {
          if (strcasecmp(config->attr_names[i], g->binattr_names[j]) == 0) {
            found = TRUE;
            config->attr_indices[i] = j;
            CONFIG_DEBUG_PRINT(("binattr %s(%s) index %u\n",
                                config->attr_param_names[i],
                                config->attr_names[i], j));
          }
        }
        if (!found) {
          fprintf(stderr, "ERROR: binary attribute %s not found\n",
                  config->attr_names[i]);
          return 1;
        }
        break;
        
      case ATTR_TYPE_CATEGORICAL:
        for (j = 0; j < g->num_catattr; j++) {
          if (strcasecmp(config->attr_names[i], g->catattr_names[j]) == 0) {
            found = TRUE;
            config->attr_indices[i] = j;
            CONFIG_DEBUG_PRINT(("catattr %s(%s) index %u\n",
                                config->attr_param_names[i],
                                config->attr_names[i], j));
          }
        }
        if (!found) {
          fprintf(stderr, "ERROR: categorical attribute %s not found\n",
                  config->attr_names[i]);
          return 1;
        }
        break;

      case ATTR_TYPE_CONTINUOUS:
        for (j = 0; j < g->num_contattr; j++) {
          if (strcasecmp(config->attr_names[i], g->contattr_names[j]) == 0) {
            found = TRUE;
            config->attr_indices[i] = j;
            CONFIG_DEBUG_PRINT(("contattr %s(%s) index %u\n",
                                config->attr_param_names[i],
                                config->attr_names[i], j));
          }
        }
        if (!found) {
          fprintf(stderr, "ERROR: continuous attribute %s not found\n",
                  config->attr_names[i]);
          return 1;
        }
        break;

      case ATTR_TYPE_SET:
        for (j = 0; j < g->num_setattr; j++) {
          if (strcasecmp(config->attr_names[i], g->setattr_names[j]) == 0) {
            found = TRUE;
            config->attr_indices[i] = j;
            CONFIG_DEBUG_PRINT(("setattr %s(%s) index %u\n",
                                config->attr_param_names[i],
                                config->attr_names[i], j));
          }
        }
        if (!found) {
          fprintf(stderr, "ERROR: set attribute %s not found\n",
                  config->attr_names[i]);
          return 1;
        }
        break;
        
      default:
        fprintf(stderr, "ERROR (internal): unknown attribute type %u\n",
                get_attr_param_type(config->attr_param_names[i]));
        return 1;
        break;
    }
  }
  return 0;
}

/*
 * build_dyadic_indices_from_names() is called after the config file
 * is parsed by parse_config_file() and also after the 
 * attributes are loaded by calling load_digraph_from_arclist_file()
 * since the relevant filenames have to be parsed from the config
 * file before their contents can be loaded. (Could enforce that
 * those filenames are earlier in the config file than the attribute
 * parameters and do it all in one pass, but did not do so).
 *
 * The attribute names specified for attribute parameter estimation
 * are searched for in the names loaded from the attribute data files
 * stored in the digraph structure g, and if they are present, the
 * indices stored in the config.dyadic_indices array. Returns nonzero on
 * error else 0.
 *
 */
int build_dyadic_indices_from_names(config_t *config,  digraph_t *g)
{
  uint_t i, j;
  bool   found;
  dyadic_type_e dyadicType;
  uint_t numGeoAttr     = 0; /* number of GeoDistance or LogGeodistance attrs */
  uint_t numEuclideanAttr = 0; /* number of EuclideanDistance attrs */
  uint_t numAttrs; /* total number of dyadic attributes */
  uint_t geoIndex              = 0;
  uint_t euclideanIndex        = 0;
  int    firstGeoIndex         = -1; /* -1 until set to first index */
  int    firstEuclideanIndex   = -1;
  char  *tmpGeoName, *tmpEuclideanName, *tmpGeoParamName, *tmpEuclideanParamName;
  dyadic_change_stats_func_t *tmpGeoFunc, *tmpEuclideanFunc;
  
  config->dyadic_indices = safe_malloc(config->num_dyadic_change_stats_funcs *
                                     sizeof(uint_t));
  config->dyadic_types = safe_malloc(config->num_dyadic_change_stats_funcs *
                                     sizeof(uint_t));  
  
  for (i = 0; i < config->num_dyadic_change_stats_funcs; i++) {
    found = FALSE;
    dyadicType = get_dyadic_param_type(config->dyadic_param_names[i]);
    switch (dyadicType) {

      case DYADIC_TYPE_GEODISTANCE:
      case DYADIC_TYPE_EUCLIDEANDISTANCE:

        for (j = 0; j < g->num_contattr; j++) {
          if (strcasecmp(config->dyadic_names[i], g->contattr_names[j]) == 0) {
            found = TRUE;
            config->dyadic_indices[i] = j;
            config->dyadic_types[i] = dyadicType;
            CONFIG_DEBUG_PRINT(("dyadic covariate type %s "
                                "contattr %s(%s) index %u\n",
                                (dyadicType == DYADIC_TYPE_GEODISTANCE ?
                                 "GEODISTANCE" :
                                 (dyadicType == DYADIC_TYPE_EUCLIDEANDISTANCE ?
                                  "EUCLIDEANDISTANCE" : "*ERROR*")),
                                config->dyadic_param_names[i],
                                config->dyadic_names[i], j));
            if (dyadicType == DYADIC_TYPE_GEODISTANCE) {
              numGeoAttr++;
            } else if (dyadicType == DYADIC_TYPE_EUCLIDEANDISTANCE) {
              numEuclideanAttr++;
            } else {
              assert(FALSE);
            }
          }
        }

        if (!found) {
          fprintf(stderr, "ERROR: dyadic covariate continuous attribute %s not found\n",
                  config->dyadic_names[i]);
          return 1;
        }
        break;

      default:
        fprintf(stderr, "ERROR (internal): unknown dyadic covariate type %u\n",
                get_dyadic_param_type(config->dyadic_param_names[i]));
        return 1;
        break;
    }
  }
  numAttrs = i;
  
  if (numAttrs > 0) {
    /* only two defined so far, GeoDistance, which requires exactly 2
       continuous attributes, for latitude and longitude respectively
       (logGeoDistance is just the same, but funciton does log of distance)
       and EuclideanDistance, which requires exactly 3 continuous attribute,
       for x, y, and z coordinates respectively.
    */
    if (numGeoAttr > 0 && numGeoAttr != 2) {
      fprintf(stderr,
              "ERROR: GeoDistance or logGeoDistance requires exactly two continuous attribute "
              "names, for latitude and longitude respectively\n");
      return 1;
    }
    if (numEuclideanAttr > 0 && numEuclideanAttr != 3) {
      fprintf(stderr,
              "ERROR: EuclideanDistance requires exactly three continuous "
              "attribute names, for x, y, and z coordinates respectively\n");
      return 1;
    }
    for (j = 0; j < numAttrs; j++) {
      switch (config->dyadic_types[j]) {
        case DYADIC_TYPE_GEODISTANCE:
          switch (geoIndex) {
            case 0:
              g->latitude_index = config->dyadic_indices[j];
              break;
            case 1:
              g->longitude_index = config->dyadic_indices[j];
              break;
            default:
              fprintf(stderr, "ERROR (internal): geoIndex == %u\n",geoIndex);
              return 1;
              break;
          }
          if (firstGeoIndex < 0) {
            firstGeoIndex = j;
          }
          geoIndex++;
          break;
          
        case DYADIC_TYPE_EUCLIDEANDISTANCE:
          switch (euclideanIndex) {
            case 0:
              g->x_index = config->dyadic_indices[j];
              break;
            case 1:
              g->y_index = config->dyadic_indices[j];
              break;
            case 2:
              g->z_index = config->dyadic_indices[j];
              break;
            default:
              fprintf(stderr, "ERROR (internal): euclideanIndex == %u\n",
                      euclideanIndex);
              return 1;
              break;
          }
          if (firstEuclideanIndex < 0) {
            firstEuclideanIndex = j;
          }
          euclideanIndex++;
          break;
          
        default:
          fprintf(stderr, "ERROR (internal): unknown dyadic type %d\n",
                  (int)config->dyadic_types[j]);
          return 1;
          break;
      }
    }
    /* Because is has two attribute names, we get two entries for the
       change statistic. But for [log]GeoDistance we only want one (it uses
       two attributes at each node), so delete the second. 
       Similarly for EuclideanDistance (But with 3 attribute names and hence
       dyadic change statistic entries).
    */
    /* TODO this system of ad hoc handling of dyadic parameters was
       alright when there was only one (or now two) types but it will
       be a real mess for any more, should make it more general (and
       hopefully simpler) */
    
    if (numGeoAttr > 0 && numEuclideanAttr > 0) {
      /* [log]GeoDistance and EuclideanDistance */
      assert(firstGeoIndex >= 0);
      assert(firstEuclideanIndex >= 0);
      assert(firstGeoIndex != firstEuclideanIndex);
      assert(get_dyadic_param_type(config->dyadic_param_names[firstGeoIndex]) == DYADIC_TYPE_GEODISTANCE);
      assert(get_dyadic_param_type(config->dyadic_param_names[firstEuclideanIndex]) == DYADIC_TYPE_EUCLIDEANDISTANCE);
      assert(config->num_dyadic_change_stats_funcs == 5);
      tmpGeoName = safe_strdup(config->dyadic_names[firstGeoIndex]);
      tmpGeoParamName = safe_strdup(config->dyadic_param_names[firstGeoIndex]);
      tmpGeoFunc = config->dyadic_change_stats_funcs[firstGeoIndex];
      tmpEuclideanName = safe_strdup(config->dyadic_names[firstEuclideanIndex]);
      tmpEuclideanParamName = safe_strdup(config->dyadic_param_names[firstEuclideanIndex]);
      tmpEuclideanFunc = config->dyadic_change_stats_funcs[firstEuclideanIndex];
      config->dyadic_names[0] = tmpGeoName;
      config->dyadic_param_names[0] = tmpGeoParamName;
      config->dyadic_change_stats_funcs[0] = tmpGeoFunc;
      config->dyadic_types[0] = DYADIC_TYPE_GEODISTANCE;
      config->dyadic_names[1] = tmpEuclideanName;
      config->dyadic_param_names[1] = tmpEuclideanParamName;
      config->dyadic_change_stats_funcs[1] = tmpEuclideanFunc;
      config->dyadic_types[1] = DYADIC_TYPE_EUCLIDEANDISTANCE;
      config->num_dyadic_change_stats_funcs = 2;      
      free(config->dyadic_names[2]);
      free(config->dyadic_names[3]);
      free(config->dyadic_names[4]);
    } else if (numGeoAttr > 0) {
      /* only [log]GeoDistance */
      assert(firstGeoIndex >= 0);
      assert(config->num_dyadic_change_stats_funcs == 2);
      assert(get_dyadic_param_type(config->dyadic_param_names[0])
             == DYADIC_TYPE_GEODISTANCE &&
             get_dyadic_param_type(config->dyadic_param_names[1])
             == DYADIC_TYPE_GEODISTANCE);
      config->num_dyadic_change_stats_funcs = 1;
      free(config->dyadic_names[1]);
    } else if (numEuclideanAttr > 0) {
      /* only EuclideanDistance */
      assert(firstEuclideanIndex >= 0);
      assert(config->num_dyadic_change_stats_funcs == 3);
      assert(get_dyadic_param_type(config->dyadic_param_names[0])
             == DYADIC_TYPE_EUCLIDEANDISTANCE &&
             get_dyadic_param_type(config->dyadic_param_names[1])
             == DYADIC_TYPE_EUCLIDEANDISTANCE &&
             get_dyadic_param_type(config->dyadic_param_names[2])
             == DYADIC_TYPE_EUCLIDEANDISTANCE);
      config->num_dyadic_change_stats_funcs = 1;
      free(config->dyadic_names[1]);
      free(config->dyadic_names[2]);
    } else {
      assert(FALSE);
    }
    
  }
  return 0;
}


/*
 * build_attr_interaction_pair_indices_from_names() is called after the
 * config file is parsed by parse_config_file() and also after the
 * attributes are loaded by calling load_digraph_from_arclist_file()
 * since the relevant filenames have to be parsed from the config file
 * before their contents can be loaded. (Could enforce that those
 * filenames are earlier in the config file than the attribute
 * parameters and do it all in one pass, but did not do so).
 *
 * The attribute names specified for attribute parameter estimation
 * are searched for in the names loaded from the attribute data files
 * stored in the digraph structure g, and if they are present, the
 * indices stored in the config.attr_indices array. Returns nonzero on
 * error else 0.
 *
 * Note that the types of attributes (binary, categorical,continuous,
 * set) are mixed together in the
 * config.attr_interaction_pairindices[] array; the index could be
 * into g->binattr or g->catattr or g->contattr or g->setattr, but
 * because each change statistics function is for one particular type,
 * it uses the index in the correct array without ambiguity.
 */
int build_attr_interaction_pair_indices_from_names(config_t *config,
                                                   const digraph_t *g)
{
  uint_t i, j, attrnum;
  bool   found;
  char  *attrname;

  config->attr_interaction_pair_indices = safe_malloc(
    config->num_attr_interaction_change_stats_funcs * sizeof(uint_pair_t));

  for (i = 0; i < config->num_attr_interaction_change_stats_funcs; i++) {
    for (attrnum = 0; attrnum < 2 /* 0=first and 1=second */; attrnum++) {
      found = FALSE;
      if (attrnum == 0) {
        attrname = config->attr_interaction_pair_names[i].first;
      } else {
        attrname = config->attr_interaction_pair_names[i].second;
      }
      switch (get_attr_interaction_param_type(
                config->attr_interaction_param_names[i])) {
        case ATTR_TYPE_BINARY:
          assert(FALSE); /* TODO binary attribute interaction */
          break;

        case ATTR_TYPE_CATEGORICAL:
          for (j = 0; j < g->num_catattr; j++) {
            if (strcasecmp(attrname, g->catattr_names[j]) == 0) {
              found = TRUE;
              if (attrnum == 0) {
                config->attr_interaction_pair_indices[i].first = j;
              } else {
                config->attr_interaction_pair_indices[i].second = j;         
              }
              CONFIG_DEBUG_PRINT(("catattr interaction [%s] %s(%s) index %u\n",
                                  attrnum == 0 ? "first" : "second",
                                  config->attr_interaction_param_names[i],
                                  attrname, j));
            }
          }
          if (!found) {
            fprintf(stderr, "ERROR: categorical attribute %s not found\n",
                    attrname);
            return 1;
          }
          break;

        case ATTR_TYPE_CONTINUOUS:
          assert(FALSE); /* TODO contiuous attribute interaction */
          break;

        case ATTR_TYPE_SET:
          assert(FALSE); /* TODO set attribute interaction (if it makes sense)*/
          break;

        default:
          fprintf(stderr, "ERROR (internal): unknown attribute type %u\n",
                  get_attr_interaction_param_type(attrname));
          return 1;
          break;
      }
    }
  }
  return 0;
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
  uint_t i;
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
  free(config->change_stats_funcs);
  free(config->param_names);
  free(config->attr_param_names);
  for (i = 0; i < config->num_attr_change_stats_funcs; i++) 
    free(config->attr_names[i]);
  free(config->attr_indices);
  for (i = 0; i < config->num_dyadic_change_stats_funcs; i++) 
    free(config->dyadic_names[i]);
  free(config->dyadic_indices);
  free(config->dyadic_types);
  free(config->attr_interaction_param_names);
  for (i = 0; i < config->num_attr_interaction_change_stats_funcs; i++)  {
    free(config->attr_interaction_pair_names[i].first);
    free(config->attr_interaction_pair_names[i].second);
  }
  free(config->attr_interaction_pair_indices);
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



/*
 * Write the allowed configuration parameters, their descriptions and 
 * default values, to stderr
 */
void dump_config_names(void)
{
  uint_t i;
  fprintf(stderr, "Configuration parameters:\n");
  for (i = 0; i < NUM_CONFIG_PARAMS; i++) {
    fprintf(stderr, "  %s: %s ", CONFIG_PARAMS[i].name,
            CONFIG_PARAMS[i].description);
    switch (CONFIG_PARAMS[i].type) {
      case PARAM_TYPE_DOUBLE:
        fprintf(stderr, "(floating point) [default %g]\n",
                *(double *)((char *)&CONFIG + CONFIG_PARAMS[i].offset));
        break;
        
      case PARAM_TYPE_UINT:
        fprintf(stderr, "(unsigned integer) [default %u]\n",
                *(uint_t *)((char *)&CONFIG + CONFIG_PARAMS[i].offset));
        break;

      case PARAM_TYPE_BOOL:
        fprintf(stderr, "(Boolean) [default %s]\n",
                *(bool *)((char *)&CONFIG + CONFIG_PARAMS[i].offset) ?
                "True" : "False");
        break;

      case PARAM_TYPE_STRING:
        fprintf(stderr, "(string)");
        if (*(char **)((char *)&CONFIG + CONFIG_PARAMS[i].offset))
          fprintf(stderr, " [default %s]\n", *(char **)((char *)&CONFIG + CONFIG_PARAMS[i].offset));
        else
          fprintf(stderr, "\n");
        break;

      case PARAM_TYPE_SET:
        fprintf(stderr, "(set of ERGM parameter names)\n");
        break;
        
      default:
      fprintf(stderr, "ERROR (internal): unknown parameter type %d\n",
              CONFIG_PARAMS[i].type);
      break;
    }
  }
}
