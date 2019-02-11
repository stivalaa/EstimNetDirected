/*****************************************************************************
 * 
 * File:    configparser.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * Parse the configuration file to get algorithm paraemeters, input filenames,
 * parameters to estimate, etc.
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
#include "configparser.h"

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
  ATTR_TYPE_CONTINUOUS      /* continuous attribute type (double) */
} attr_type_e;

/* ERGM dyadic covariate parameter type */
typedef enum dyadic_type_e {
  DYADIC_TYPE_INVALID,       /* invalid type, used as error return value */
  DYADIC_TYPE_GEODISTANCE    /* continuous geographic distance from lat/long */
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

/*****************************************************************************
 *
 * constants
 *
 ****************************************************************************/

static const size_t BUFSIZE = 16384;  /* line buffer size for reading files */
static const size_t TOKSIZE = 1024;   /* maximum size of a token */
static const char   COMMENT_CHAR = '#'; /* comment character */
static const char   OPEN_SET_CHAR = '{';  /* set of parameter vals open */
static const char   CLOSE_SET_CHAR = '}'; /* set of parameter vals close */
static const char   OPEN_PAREN_CHAR = '(';
static const char   CLOSE_PAREN_CHAR = ')';


/* True and False values for Boolean config value type. Not case sensitive */
const char *TRUE_STR = "true";
const char *FALSE_STR = "false";

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
  "binary/categorical/continuous attribute parameters to estimate"},

  {DYADIC_PARAMS_STR, PARAM_TYPE_SET,       0, /*no offset, coded explicitly*/
  "dyadic covariate parameters to estimate"},

  {ATTR_INTERACTION_PARAMS_STR,PARAM_TYPE_SET, 0,/*no offset, coded explicitly*/
   "attribute pair interaction parameters to estimate"}
};
static const uint_t NUM_CONFIG_PARAMS = sizeof(CONFIG_PARAMS) /
  sizeof(CONFIG_PARAMS[0]);

/*
 * Structural parameters allowed as the names in the set for the 
 * structParams parameter. Names are not case sensitive 
 */
static const struct_param_t STRUCT_PARAMS[] =
{
  {ARC_PARAM_STR,       changeArc},
  {"Reciprocity",       changeReciprocity},
  {"Sink",              changeSink},
  {"Source",            changeSource},
  {"Isolates",          changeIsolates},
  {"TwoPaths",          changeTwoPath},
  {"InTwoStars",        changeInTwoStars},
  {"OutTwoStars",       changeOutTwoStars},
  {"TransitiveTriangles",changeTransitiveTriad},
  {"CyclicTriangles",   changeCyclicTriad},
  {"AltInStars",        changeAltInStars},
  {"AltOutStars",       changeAltOutStars},
  {"AltKTrianglesT",    changeAltKTrianglesT},
  {"AltKTrianglesC",    changeAltKTrianglesC},
  {"AltKTrianglesD",    changeAltKTrianglesD},
  {"AltKTrianglesU",    changeAltKTrianglesU},  
  {"AltTwoPathsT",      changeAltTwoPathsT},
  {"AltTwoPathsD",      changeAltTwoPathsD},
  {"AltTwoPathsU",      changeAltTwoPathsU},
  {"AltTwoPathsTD",     changeAltTwoPathsTD}
};
static const uint_t NUM_STRUCT_PARAMS = sizeof(STRUCT_PARAMS) /
  sizeof(STRUCT_PARAMS[0]);

/*
 * Attribute parameters allowed as the names in the set for the 
 * attrParams  parameters. Names are not case sensitive.
 */
static const attr_param_t ATTR_PARAMS[] =
{
  {"Sender",       ATTR_TYPE_BINARY,      changeSender},
  {"Receiver",     ATTR_TYPE_BINARY,      changeReceiver},
  {"Interaction",  ATTR_TYPE_BINARY,      changeInteraction},
  {"Matching",     ATTR_TYPE_CATEGORICAL, changeMatching},
  {"MatchingReciprocity",    ATTR_TYPE_CATEGORICAL, changeMatchingReciprocity},
  {"Mismatching",            ATTR_TYPE_CATEGORICAL, changeMismatching},
  {"MismatchingReciprocity", ATTR_TYPE_CATEGORICAL, changeMismatchingReciprocity},
  {"ContinuousSender",       ATTR_TYPE_CONTINUOUS, changeContinuousSender},
  {"ContinuousReceiver",     ATTR_TYPE_CONTINUOUS, changeContinuousReceiver},
  {"Diff",                   ATTR_TYPE_CONTINUOUS, changeDiff},
  {"DiffReciprocity",        ATTR_TYPE_CONTINUOUS, changeDiffReciprocity},
  {"DiffSign",               ATTR_TYPE_CONTINUOUS, changeDiffSign}
};
static const uint_t NUM_ATTR_PARAMS = sizeof(ATTR_PARAMS) /
  sizeof(ATTR_PARAMS[0]);


/*
 * Dyadic covariate parameters allowed sa the names in the set for the
 * dyadicParams parameters. Names are not case sensitive.
 */
static const dyadic_param_t DYADIC_PARAMS[] =
{
  {"GeoDistance",    DYADIC_TYPE_GEODISTANCE,   changeGeoDistance},
  {"logGeoDistance", DYADIC_TYPE_GEODISTANCE,   changeLogGeoDistance}
};
static const uint_t NUM_DYADIC_PARAMS = sizeof(DYADIC_PARAMS) /
  sizeof(DYADIC_PARAMS[0]);

/*
 * Attribute pair interaction parameters allowed as the names in the
 * set for the attrInteractionParams parameters. Names are not case sensitive.
 */
static const attr_interaction_param_t ATTR_INTERACTION_PARAMS[] =
{
  {"MatchingInteraction",     ATTR_TYPE_CATEGORICAL, changeMatchingInteraction},
};
static const uint_t NUM_ATTR_INTERACTION_PARAMS =
  sizeof(ATTR_INTERACTION_PARAMS) / sizeof(ATTR_INTERACTION_PARAMS[0]);


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
 * Skip whitespace in (open read) file f. First non-whitespace char
 * is returned (or EOF on end or error).
 */
static int fskip(FILE *f)
{
  int c = fgetc(f);
  while (c != EOF && isspace(c))
    c = fgetc(f);
  return c;
}

/* Return nonzero if c is a standalone token character else 0 */
static int isSingleCharToken(int c)
{
  return c == '=' || c == ',' ||
    c == OPEN_PAREN_CHAR || c == CLOSE_PAREN_CHAR ||
    c == OPEN_SET_CHAR || c == CLOSE_SET_CHAR;
}

/* Return nonzero if c is a character allowed in a token,
   but not a token on its own, else 0 */
static int istokenchar(int c)
{
  return c != COMMENT_CHAR && (isalnum(c) || ispunct(c)) &&
    !isSingleCharToken(c);
}

/* return nonzero if c is a char allowed in a parameter name else 0 */
static int isParamNameChar(int c)
{
  return isalnum(c) || c == '_';
}

/* return nonzero if s is a lexically valid config parameter name else 0 */
static int isValidParamname(const char *s)
{
  const char *p;
  for (p = s; p != '\0'; p++)
    if (!isParamNameChar(*p))
      return 1;
  return 0;
}

/*
 * get next token in infile (open for read). Returns NULL if none. The
 * token is put in the token parameter (which is return value). token
 * must be allocated by caller with space for at least TOKSIZE chars.
 */
static char *get_token(FILE *infile, char *token)
{
  char buf[BUFSIZE];
  int   c; /* unsigned char cast to int from getc() etc. */
  int   i = 0;
  c = fskip(infile);
  while (c != EOF && c == COMMENT_CHAR) {
    /* keep disarding lines starting with comment char until we get a token*/
    fgets(buf, sizeof(buf), infile);
    c = fskip(infile);
  }
  if (c == EOF)
    return NULL;
  if (isSingleCharToken(c)) { /* '=' is always a token all on its own */
    token[i++] = c;
    token[i] = '\0';
    return token;
  }
  while (c != EOF && istokenchar(c) && i < (int)TOKSIZE-1) {
    token[i++] = c;
    c = fgetc(infile);
  }
  token[i] = '\0';
  ungetc(c, infile);
  return token;
}


/* 
 * Get simple "paramname = value" pair from the infile
 * and return them in the output parameters paramname and value,
 * which must be allocated by caller with space for at least TOKSIZE chars
 * each. Return 0 if succesfully got parameter and value, 
 * +ve if did not (but no error), -ve on error.
 */
static int get_paramname_value(FILE *infile, char *paramname, char *value)
{
  char      tokenbuf[TOKSIZE];
  char     *token;
  
  if (!(token = get_token(infile, tokenbuf))) {
    return 1; /* not necessarily an error, may be no tokens left in file */
  }
  if (!isValidParamname(token)) {
    fprintf(stderr, "ERROR: invalid config parameter '%s'\n", token);
    return -1;
  }
  strncpy(paramname, token, TOKSIZE);
  if (!(token = get_token(infile, tokenbuf)) || strcmp(token, "=") != 0) {
    fprintf(stderr, "ERROR: expecting '=' after configuration parameter "
            "%s but found '%s'\n", paramname, token);
    return -1;
  }
  if (!(token = get_token(infile, tokenbuf))) {
    fprintf(stderr, "ERROR: Did not find value for configuration parameter %s\n",
            paramname);
    return -1;
  }
  strncpy(value, token, TOKSIZE);
  return 0;
}


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
        CONFIG_DEBUG_PRINT(("attrInteractionParam %s('%s')\n",
                            paramName, token));
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
          CONFIG.num_attr_interaction_change_stats_funcs++;
          num_attr_names++;
        }  else if (num_attr_names == 1) {
          CONFIG.attr_interaction_pair_names[
            CONFIG.num_attr_interaction_change_stats_funcs].second = safe_strdup(token);
          num_attr_names++;
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
  free(config->attr_interaction_param_names);
  for (i = 0; i < config->num_attr_interaction_change_stats_funcs; i++)  {
    free(config->attr_interaction_pair_names[i].first);
    free(config->attr_interaction_pair_names[i].second);
  }
  free(config->attr_interaction_pair_indices);
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

/*
 * Given name of attribute parameter, returns its type, by searching
 * in the ATTR_PARAMS const array. Note this is a linear search but
 * does not matter as a small hand-coded array anyway.
 */
static attr_type_e get_attr_param_type(const char *attrParamName)
{
  uint_t i;

  for (i = 0; i < NUM_ATTR_PARAMS; i++) {
    if (strcasecmp(ATTR_PARAMS[i].name, attrParamName) == 0) {
      return ATTR_PARAMS[i].type;
    }
  }
  return ATTR_TYPE_INVALID;
}

/*
 * Given name of dyadic covariate parameter, returns its type, by searching
 * in the DYADIC_PARAMS const array. Note this is a linear search but
 * does not matter as a small hand-coded array anyway.
 */
static dyadic_type_e get_dyadic_param_type(const char *dyadicParamName)
{
  uint_t i;

  for (i = 0; i < NUM_DYADIC_PARAMS; i++) {
    if (strcasecmp(DYADIC_PARAMS[i].name, dyadicParamName) == 0) {
      return DYADIC_PARAMS[i].type;
    }
  }
  return DYADIC_TYPE_INVALID;
}

/*
 * Given name of attribute interaction parameter, returns its type, by searching
 * in the ATTR_INTERACTION_PARAMS const array. Note this is a linear search but
 * does not matter as a small hand-coded array anyway.
 */
static attr_type_e get_attr_interaction_param_type(const char
                                                   *attrInteractionParamName)
{
  uint_t i;

  for (i = 0; i < NUM_ATTR_INTERACTION_PARAMS; i++) {
    if (strcasecmp(ATTR_INTERACTION_PARAMS[i].name,
                   attrInteractionParamName) == 0) {
      return ATTR_INTERACTION_PARAMS[i].type;
    }
  }
  return ATTR_TYPE_INVALID;
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
 * Note that the types of attributes (binary, categorical,continuous)
 * are mixed together in the config.attr_indices[] array; the index
 * could be into g->binattr or g->catattr or g->contattr, but because
 * each change statistics function is for one particular type, it uses
 * the index in the correct array without ambiguity.
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
 * Note that currently only the DYADIC_TYPE_GEODISTANCE is used
 * and GeoDistance (or logGeoDistance) is the only dyadic covariate parameter: this
 * is a special case which takes exactly two continuous attribute
 * names and sets their indices as the latitude_index and 
 * longitude_index in the digraph structure.
 *
 */
int build_dyadic_indices_from_names(config_t *config,  digraph_t *g)
{
  uint_t i, j;
  bool   found;

  config->dyadic_indices = safe_malloc(config->num_dyadic_change_stats_funcs *
                                     sizeof(uint_t));
  
  for (i = 0; i < config->num_dyadic_change_stats_funcs; i++) {
    found = FALSE;
    switch (get_dyadic_param_type(config->dyadic_param_names[i])) {

      case DYADIC_TYPE_GEODISTANCE:
        for (j = 0; j < g->num_contattr; j++) {
          if (strcasecmp(config->dyadic_names[i], g->contattr_names[j]) == 0) {
            found = TRUE;
            config->dyadic_indices[i] = j;
            CONFIG_DEBUG_PRINT(("dyadic covariate contattr %s(%s) index %u\n",
                                config->dyadic_param_names[i],
                                config->dyadic_names[i], j));
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

  if (i > 0) {
    /* only one defined so far, GeoDistance, which requires exactly 2
       continuous attributes, for latitude and longitude respectively
       (logGeoDistance is just the same, but funciton does log of distance)*/
    if (i != 2) {
      fprintf(stderr,
              "ERROR: GeoDistance or logGeoDistance requires exactly two continuous attribute "
              "names, for latitude and longitude respectively\n");
      return 1;
    }
    g->latitude_index = config->dyadic_indices[0];
    g->longitude_index = config->dyadic_indices[1];
    /* Because is has two attribute names, we get two entries for the
       change statistic. But for [log]GeoDistance we only want one (it uses
       two attributes at each node), so delete the second. */
    assert(config->num_dyadic_change_stats_funcs == 2);
    assert(get_dyadic_param_type(config->dyadic_param_names[0])
           == DYADIC_TYPE_GEODISTANCE &&
      get_dyadic_param_type(config->dyadic_param_names[1])
           == DYADIC_TYPE_GEODISTANCE);
    config->num_dyadic_change_stats_funcs = 1;
    free(config->dyadic_names[1]);
  }
  return 0;
}



/*
 * Write the allowed ERGM estimation parameters to stderr 
 */
void dump_parameter_names(void)
{
  uint_t i;
  fprintf(stderr, "Structural parameters (%s):\n", STRUCT_PARAMS_STR);
  for (i = 0; i < NUM_STRUCT_PARAMS; i++) {
    fprintf(stderr, "  %s\n", STRUCT_PARAMS[i].name);
  }
  fprintf(stderr, "Attribute parameters (%s):\n", ATTR_PARAMS_STR);
  for (i = 0; i < NUM_ATTR_PARAMS; i++) {
    fprintf(stderr, "  %s (%s)\n", ATTR_PARAMS[i].name,
            ATTR_PARAMS[i].type == ATTR_TYPE_BINARY ? "binary" :
            (ATTR_PARAMS[i].type == ATTR_TYPE_CATEGORICAL ? "categorical" :
             (ATTR_PARAMS[i].type == ATTR_TYPE_CONTINUOUS ? "continuous" :
              "*UNKNOWN*")));
  }
  fprintf(stderr, "Dyadic covariate parameters (%s):\n", DYADIC_PARAMS_STR);
  for (i = 0; i < NUM_DYADIC_PARAMS; i++) {
    fprintf(stderr, " %s (%s)\n", DYADIC_PARAMS[i].name,
            DYADIC_PARAMS[i].type == DYADIC_TYPE_GEODISTANCE ?
            "latitude,longitude" : "*UNKNOWN*");
  }
  fprintf(stderr, "Attribute interaction parameters (%s)\n",
          ATTR_INTERACTION_PARAMS_STR);
  for (i = 0; i < NUM_ATTR_INTERACTION_PARAMS; i++) {
    fprintf(stderr, " %s (%s)\n", ATTR_INTERACTION_PARAMS[i].name,
            ATTR_INTERACTION_PARAMS[i].type == ATTR_TYPE_BINARY ? "binary" :
            (ATTR_INTERACTION_PARAMS[i].type == ATTR_TYPE_CATEGORICAL ? "categorical" :
             (ATTR_INTERACTION_PARAMS[i].type == ATTR_TYPE_CONTINUOUS ? "continuous" :
              "*UNKNOWN*")));
  }
  fprintf(stderr, "\n");
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
