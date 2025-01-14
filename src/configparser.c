/*****************************************************************************
 * 
 * File:    configparser.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * Parse the configuration file to get algorithm paraemeters, input filenames,
 * parameters to estimate, etc.
 *
 * This file for constants, functions etc. shared by estimation and simulation
 * (and potentially other) configuration parsing.
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
 * externally visible constant definitions
 *
 ****************************************************************************/

const size_t TOKSIZE = 8192;   /* maximum size of a token */

/*****************************************************************************
 *
 * local constant definitions
 *
 ****************************************************************************/

static const size_t BUFSIZE = 16384;  /* line buffer size for reading files */
static const char   COMMENT_CHAR = '#'; /* comment character */
static const char   OPEN_SET_CHAR = '{';  /* set of parameter vals open */
static const char   CLOSE_SET_CHAR = '}'; /* set of parameter vals close */
static const char   OPEN_PAREN_CHAR = '(';
static const char   CLOSE_PAREN_CHAR = ')';
static const char   QUOTE_STR_CHAR   = '"'; /* string literal open/close */

/* True and False values for Boolean config value type. Not case sensitive */
static const char *TRUE_STR = "true";
static const char *FALSE_STR = "false";


/*** wide tables, edit in 132 column terminal window  ***/

/*
 * Structural parameters allowed as the names in the set for the 
 * structParams parameter. Names are not case sensitive 
 */
static const struct_param_t STRUCT_PARAMS[] =
{
  {ARC_PARAM_STR,          STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeArc},
  {"Reciprocity",          STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeReciprocity},
  {"Sink",                 STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeSink},
  {"Source",               STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeSource},
  {"Isolates",             STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_BOTH,       NETWORK_MODE_BOTH,    changeIsolates},
  {"TwoPaths",             STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_BOTH,       NETWORK_MODE_BOTH,    changeTwoPath},
  {"InTwoStars",           STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeInTwoStars},
  {"OutTwoStars",          STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeOutTwoStars},
  {"TransitiveTriangles",  STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeTransitiveTriad},
  {"CyclicTriangles",      STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeCyclicTriad},
  {"AltInStars",           STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeAltInStars},
  {"AltOutStars",          STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeAltOutStars},
  {"AltKTrianglesT",       STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeAltKTrianglesT},
  {"AltKTrianglesC",       STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeAltKTrianglesC},
  {"AltKTrianglesD",       STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeAltKTrianglesD},
  {"AltKTrianglesU",       STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeAltKTrianglesU},
  {"AltTwoPathsT",         STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeAltTwoPathsT},
  {"AltTwoPathsD",         STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeAltTwoPathsD},
  {"AltTwoPathsU",         STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeAltTwoPathsU},
  {"AltTwoPathsTD",        STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeAltTwoPathsTD},
  {"Loop",                 STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_BOTH,       NETWORK_MODE_ONEMODE, changeLoop},
  {"LoopInteraction",      STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeLoopInteraction},
  {EDGE_PARAM_STR,         STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_BOTH,    changeEdge},
  {"TwoStars",             STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_ONEMODE, changeTwoStars},
  {"AltStars",             STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_ONEMODE, changeAltStars},
  {"AltTwoPaths",          STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_BOTH,    changeAltTwoPaths},
  {"AltKTriangles",        STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_ONEMODE, changeAltKTriangles},
  {"BipartiteTwoStarsA",   STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteTwoStarsA},
  {"BipartiteTwoStarsB",   STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteTwoStarsB},
  {"BipartiteThreeStarsA", STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteThreeStarsA},
  {"BipartiteThreeStarsB", STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteThreeStarsB},
  {"BipartiteAltStarsA",   STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteAltStarsA},
  {"BipartiteAltStarsB",   STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteAltStarsB},
  {"FourCycles"        ,   STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_BOTH,    changeFourCycles},
  {"ThreePaths",           STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_BOTH,    changeThreePaths},
  {"IsolateEdges",         STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_BOTH,    changeIsolateEdges},
  {"BipartiteAltKCyclesA", STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteAltKCyclesA},
  {"BipartiteAltKCyclesB", STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteAltKCyclesB},
{"BipartiteAltK4CyclesA",STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteAltK4CyclesA},
{"BipartiteAltK4CyclesB",STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteAltK4CyclesB},
  {"BipartiteIsolatesA",   STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteIsolatesA},
  {"BipartiteIsolatesB",   STRUCT_PARAM_TYPE_NONE,   NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteIsolatesB},
  {"FourCyclesNodePower",  STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_BOTH,    changePowerFourCycles},
  {"BipartiteFourCyclesNodePowerA", STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartitePowerFourCyclesA},
  {"BipartiteFourCyclesNodePowerB", STRUCT_PARAM_TYPE_LAMBDA, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartitePowerFourCyclesB}
};
static const uint_t NUM_STRUCT_PARAMS = sizeof(STRUCT_PARAMS) /
  sizeof(STRUCT_PARAMS[0]);


/*
 * Attribute parameters allowed as the names in the set for the 
 * attrParams  parameters. Names are not case sensitive.
 */
static const attr_param_t ATTR_PARAMS[] =
{
  {"Sender",                      ATTR_TYPE_BINARY,      ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeSender},
  {"Receiver",                    ATTR_TYPE_BINARY,      ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeReceiver},
  {"Interaction",                 ATTR_TYPE_BINARY,      ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_BOTH,       NETWORK_MODE_BOTH,    changeInteraction},
  {"Matching",                    ATTR_TYPE_CATEGORICAL, ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_BOTH,       NETWORK_MODE_BOTH,    changeMatching},
  {"MatchingReciprocity",         ATTR_TYPE_CATEGORICAL, ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeMatchingReciprocity},
  {"Mismatching",                 ATTR_TYPE_CATEGORICAL, ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_BOTH,       NETWORK_MODE_BOTH,    changeMismatching},
  {"MismatchingReciprocity",      ATTR_TYPE_CATEGORICAL, ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeMismatchingReciprocity},
  {"MismatchingTransitiveTriad",  ATTR_TYPE_CATEGORICAL, ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeMismatchingTransitiveTriad},
  {"MismatchingTransitiveTies",   ATTR_TYPE_CATEGORICAL, ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeMismatchingTransitiveTies},
  {"ContinuousSender",            ATTR_TYPE_CONTINUOUS,  ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeContinuousSender},
  {"ContinuousReceiver",          ATTR_TYPE_CONTINUOUS,  ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeContinuousReceiver},
  {"Diff",                        ATTR_TYPE_CONTINUOUS,  ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_BOTH,       NETWORK_MODE_BOTH,    changeDiff},
  {"Sum",                         ATTR_TYPE_CONTINUOUS,  ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_BOTH,       NETWORK_MODE_BOTH,    changeSum},
  {"DiffReciprocity",             ATTR_TYPE_CONTINUOUS,  ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeDiffReciprocity},
  {"DiffSign",                    ATTR_TYPE_CONTINUOUS,  ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeDiffSign},
  {"DiffDirSR",                   ATTR_TYPE_CONTINUOUS,  ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeDiffDirSR},
  {"DiffDirRS",                   ATTR_TYPE_CONTINUOUS,  ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_DIRECTED,   NETWORK_MODE_ONEMODE, changeDiffDirRS},
  {"JaccardSimilarity",           ATTR_TYPE_SET,         ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_BOTH,       NETWORK_MODE_BOTH,    changeJaccardSimilarity},
  {"Activity",                    ATTR_TYPE_BINARY,      ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_BOTH,    changeActivity},
  {"BipartiteActivityA",          ATTR_TYPE_BINARY,      ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteActivityA},
  {"BipartiteActivityB",          ATTR_TYPE_BINARY,      ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteActivityB},
  {"BipartiteContinuousActivityA",ATTR_TYPE_CONTINUOUS,  ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteContinuousActivityA},
  {"BipartiteContinuousActivityB",ATTR_TYPE_CONTINUOUS,  ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteContinuousActivityB},
  {"BipartiteTwoPathSumA",        ATTR_TYPE_CONTINUOUS,  ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteTwoPathSumA},
  {"BipartiteTwoPathSumB",        ATTR_TYPE_CONTINUOUS,  ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteTwoPathSumB},
  {"BipartiteTwoPathDiffA",       ATTR_TYPE_CONTINUOUS,  ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteTwoPathDiffA},
  {"BipartiteTwoPathDiffB",       ATTR_TYPE_CONTINUOUS,  ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteTwoPathDiffB},
  {"BipartiteTwoPathMatchingA",   ATTR_TYPE_CATEGORICAL, ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteTwoPathMatchingA},
  {"BipartiteTwoPathMatchingB",   ATTR_TYPE_CATEGORICAL, ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteTwoPathMatchingB},
  {"BipartiteTwoPathMismatchingA",ATTR_TYPE_CATEGORICAL, ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteTwoPathMatchingA},
  {"BipartiteTwoPathMismatchingB",ATTR_TYPE_CATEGORICAL, ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteTwoPathMatchingB},
  {"BipartiteNodematchAlphaA",    ATTR_TYPE_CATEGORICAL, ATTR_EXP_TYPE_EXPONENT, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteNodematchAlphaA},
  {"BipartiteNodematchBetaA",     ATTR_TYPE_CATEGORICAL, ATTR_EXP_TYPE_EXPONENT, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteNodematchBetaA},
  {"BipartiteNodematchAlphaB",    ATTR_TYPE_CATEGORICAL, ATTR_EXP_TYPE_EXPONENT, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteNodematchAlphaB},
  {"BipartiteNodematchBetaB",     ATTR_TYPE_CATEGORICAL, ATTR_EXP_TYPE_EXPONENT, NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteNodematchBetaB},
  {"BipartiteExactlyOneNeighbourA",ATTR_TYPE_BINARY,     ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteExactlyOneNeighbourA},
  {"BipartiteExactlyOneNeighbourB",ATTR_TYPE_BINARY,     ATTR_EXP_TYPE_NONE,     NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_TWOMODE, changeBipartiteExactlyOneNeighbourB}
};
static const uint_t NUM_ATTR_PARAMS = sizeof(ATTR_PARAMS) /
  sizeof(ATTR_PARAMS[0]);


/*
 * Dyadic covariate parameters allowed as the names in the set for the
 * dyadicParams parameters. Names are not case sensitive.
 */
static const dyadic_param_t DYADIC_PARAMS[] =
{
  {"GeoDistance",       DYADIC_TYPE_GEODISTANCE,       NETWORK_TYPE_BOTH, NETWORK_MODE_BOTH, changeGeoDistance},
  {"logGeoDistance",    DYADIC_TYPE_GEODISTANCE,       NETWORK_TYPE_BOTH, NETWORK_MODE_BOTH, changeLogGeoDistance},
  {"EuclideanDistance", DYADIC_TYPE_EUCLIDEANDISTANCE, NETWORK_TYPE_BOTH, NETWORK_MODE_BOTH, changeEuclideanDistance}
};
static const uint_t NUM_DYADIC_PARAMS = sizeof(DYADIC_PARAMS) /
  sizeof(DYADIC_PARAMS[0]);

/*
 * Attribute pair interaction parameters allowed as the names in the
 * set for the attrInteractionParams parameters. Names are not case sensitive.
 */
static const attr_interaction_param_t ATTR_INTERACTION_PARAMS[] =
{
  {"MatchingInteraction",     ATTR_TYPE_CATEGORICAL, NETWORK_TYPE_BOTH,       NETWORK_MODE_BOTH, changeMatchingInteraction},
  {"BinaryPairInteraction",   ATTR_TYPE_BINARY,      NETWORK_TYPE_UNDIRECTED, NETWORK_MODE_BOTH, changeBinaryPairInteraction}
};
static const uint_t NUM_ATTR_INTERACTION_PARAMS =
  sizeof(ATTR_INTERACTION_PARAMS) / sizeof(ATTR_INTERACTION_PARAMS[0]);



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
 * If struct_param_type for the name parsed is
 *  STRUCT_PARAM_TYPE_LAMBDA then the parameter
 * can optionally take a lambda (decay parameter) value in parentheses e.g.
 * AltKTrianglesT(2.0). This lambda values are set in the
 * config struct in the param_lambda list. If not specified then set to
 * the default value. For these with no lambda it is just set to 0, an
 * invalid value,  and not used (but placeholder so lists line up).
 * If requireErgmValue is TRUE then ERGM parameters require a value
 * (supplied with name = value format) for use in simulation (otherwise
 * no value allowed, used for estimation).
 * Return nonzero on error else zero.
 */
static int parse_struct_params(FILE *infile, param_config_t *pconfig,
                               bool requireErgmValue)
{
  char        tokenbuf[TOKSIZE];
  char       *token;
  bool        found_paramname;
  bool        last_token_was_paramname = FALSE;
  uint_t      i;
  char        *endptr; /* for strtod() */
  char        paramname[TOKSIZE];  /* parameter name buffer */
  double      value = 0;
  double      lambda_value = 0; /* valid values are > 1.0 */
  bool        got_token_after_paramname = FALSE;
  
  if (!(token = get_token(infile, tokenbuf))) {
    fprintf(stderr, "ERROR: no tokens for structParams\n");
    return 1; 
  }
  CONFIG_DEBUG_PRINT(("parse_struct_params [1] token '%s'\n", token));  
  while (token && !(strlen(token) == 1 && token[0] == CLOSE_SET_CHAR)) {
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
      strncpy(paramname, token, TOKSIZE);

      if (STRUCT_PARAMS[i].struct_param_type == STRUCT_PARAM_TYPE_LAMBDA) {
        /* The paramter can optionally take a decay (lambda) value in
           parentheses, e.g.  AltKTrianglesT(2.0) */
        CONFIG_DEBUG_PRINT(("  [may take lambda]\n"));
        if (!(token = get_token(infile, tokenbuf))) {
          fprintf(stderr, "ERROR: no tokens after structParam %s\n", paramname);
          return 1; 
        }
        CONFIG_DEBUG_PRINT(("parse_struct_params [2] token '%s'\n", token));
        got_token_after_paramname = TRUE;
        if (strlen(token) == 1 && token[0] == OPEN_PAREN_CHAR) {
          if (!(token = get_token(infile, tokenbuf))) {
            fprintf(stderr, "ERROR: expecting %c <lambda value> %c "
                    "after structParam %s but found no token'\n",
                    OPEN_PAREN_CHAR, CLOSE_PAREN_CHAR, paramname);
            return 1;
          }
          CONFIG_DEBUG_PRINT(("parse_struct_params [3] token '%s'\n", token));
          got_token_after_paramname = FALSE;
          lambda_value = strtod(token, &endptr);
          if (*endptr != '\0') {
            fprintf(stderr, "ERROR: expecting floating point"
                    " value for <lambda value> in %c <lambda value> %c "
                    "after structParam %s but found '%s'\n",
                    OPEN_PAREN_CHAR, CLOSE_PAREN_CHAR, paramname, token);
            return 1;
          }
          if (lambda_value <= 1.0) {
            fprintf(stderr, "ERROR: lambda value must be > 1.0 but got value "
                            "%g for structParam %s\n", lambda_value, paramname);
            return 1;
          }
          if (!(token = get_token(infile, tokenbuf))) {
            fprintf(stderr, "ERROR: expecting '%c' after lambda value %g for"
                    " structParam %s, but on tokens found\n",
                    CLOSE_PAREN_CHAR, lambda_value, paramname);
            return 1; 
          }
          CONFIG_DEBUG_PRINT(("parse_struct_params [4] token '%s'\n", token));
          CONFIG_DEBUG_PRINT(("%s lambda = %g\n", paramname, lambda_value));
        } else {
          /* no open paren char found so lambda not specified */
          lambda_value = DEFAULT_LAMBDA; /* so use default value */
          CONFIG_DEBUG_PRINT(("%s default lambda = %g\n", paramname,
                              lambda_value));
          assert(lambda_value > 1.0); /* default value must be valid */
        }
      } else {
        lambda_value = 0; /* use 0 as placeholder if no lambda for this param*/
      }

      if (requireErgmValue) {
        if (!got_token_after_paramname) {
          if (!(token = get_token(infile, tokenbuf))) {
            fprintf(stderr, "ERROR: structParams expecting 'name = value' pairs separated by comma (%s)\n", paramname);
            return -1;
          }
          CONFIG_DEBUG_PRINT(("parse_struct_params [5] token '%s'\n", token));
        }
        got_token_after_paramname = FALSE;
        if (strcmp(token, "=") != 0) {
          fprintf(stderr, "ERROR: structParams expecting 'name = value' pairs separated by comma (%s)\n", paramname);
          return 1;
        }
        if (!(token = get_token(infile, tokenbuf))) {
          fprintf(stderr, "ERROR: Did not find value for structParams %s\n",
                  paramname);
          return -1;
        }
        CONFIG_DEBUG_PRINT(("parse_struct_params [6] token '%s'\n", token));        
        value = strtod(token, &endptr);
        if (*endptr != '\0') {
          fprintf(stderr, "ERROR: expecting floating point value for structParam %s but got '%s'\n", paramname, token);
          return 1;
        }
        CONFIG_DEBUG_PRINT(("structParam value %g\n", value));
      }
      
      pconfig->param_names = (const char **)safe_realloc(pconfig->param_names,
                                        (pconfig->num_change_stats_funcs + 1) *
                                        sizeof(const char *));
      pconfig->change_stats_funcs = (change_stats_func_t **)
        safe_realloc(pconfig->change_stats_funcs,
                  (pconfig->num_change_stats_funcs + 1) *
                     sizeof(change_stats_func_t *));
      pconfig->param_lambdas = (double *)safe_realloc(pconfig->param_lambdas,
                       (pconfig->num_change_stats_funcs + 1) * sizeof(double));
      pconfig->param_names[pconfig->num_change_stats_funcs] =
        STRUCT_PARAMS[i].name;
      pconfig->change_stats_funcs[pconfig->num_change_stats_funcs] =
        STRUCT_PARAMS[i].change_stats_func;
      pconfig->param_lambdas[pconfig->num_change_stats_funcs] = lambda_value;
      if (requireErgmValue) {
        pconfig->param_values = (double *)safe_realloc(pconfig->param_values,
                         (pconfig->num_change_stats_funcs + 1) * sizeof(double));
        pconfig->param_values[pconfig->num_change_stats_funcs] = value;
      }
      pconfig->num_change_stats_funcs++;
    }
    CONFIG_DEBUG_PRINT(("got_token_after_paramname = %d\n",
                        got_token_after_paramname));
    if (!got_token_after_paramname) {
      token = get_token(infile, tokenbuf);
      CONFIG_DEBUG_PRINT(("parse_struct_params [7] token '%s'\n", token));
    }
    got_token_after_paramname = FALSE;
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
 *
 * If requireErgmValue is TRUE then ERGM parameters require a value
 * (supplied with name = value format) for use in simulation (otherwise
 * no value allowed, used for estimation).
 *
 * If the attr_exponent_type is ATTR_EXP_TYPE_EXPONENT then the
 * attribute name must be followed by an exponent value (a floating
 * point value in the interval [0, 1]) as used by some change
 * stasistics, specifically for example BipartiteNodematchAlphaA and
 * BipartiteNodematchBetaA for the alpha and beta values
 * respectively. e.g.: "BipartiteNodematchBetaA(gender, 0.1)"
 */
static int parse_one_attr_param(const char *paramName,
                                attr_change_stats_func_t *attr_change_stats_func,
                                FILE *infile, param_config_t *pconfig,
                                bool requireErgmValue,
                                attr_exponent_type_e exponent_type)
{
  char      tokenbuf[TOKSIZE];
  char     *token;
  bool      last_token_was_attrname = FALSE;
  bool      opening = TRUE; /* true for first iteration only to expect '(' */
  double    value   = 0;
  char     *endptr; /* for strtod() */
  char      attrname[TOKSIZE];  /* attribute name buffer */
  double    exponent_value = -1; /* valid values are >= 0 */
  
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
        strncpy(attrname, token, TOKSIZE);

        if (exponent_type == ATTR_EXP_TYPE_EXPONENT) {
          /* The attribute name must be followed by floating
             point value in interval [0, 1] used e.g. as an exponent
             or decay parameter */
          CONFIG_DEBUG_PRINT(("attrParam %s('%s') expecting exponent\n",
                              paramName, attrname));
          if (!(token = get_token(infile, tokenbuf)) ||
              strcmp(token, ",") != 0) {
            fprintf(stderr, "attrParam %s('%s') expecting comma then exponent\n",
                    paramName, attrname);
            return 1;
          }
          if (!(token = get_token(infile, tokenbuf))) {
            fprintf(stderr, "attrParam %s('%s') expecting exponent\n",
                    paramName, attrname);
            return 1;
          }
          CONFIG_DEBUG_PRINT(("parse_one_attr_param token '%s'\n", token));        
          exponent_value = strtod(token, &endptr);
          if (*endptr != '\0') {
            fprintf(stderr, "ERROR: expecting floating point value for attrParam %s(%s) exponent but got '%s'\n", paramName, attrname, token);
            return 1;
          }
          CONFIG_DEBUG_PRINT(("attrParam %s(%s) exponent %g\n", paramName, attrname, exponent_value));
          if (exponent_value < 0.0 || exponent_value > 1.0) {
            fprintf(stderr, "ERROR: attrParam %s(%s) exponent value %g is not in interval [0, 1]\n",
                    paramName, attrname, exponent_value);
            return 1;
          }
        }
                

        if (requireErgmValue) {
          if (!(token = get_token(infile, tokenbuf))) {
            fprintf(stderr, "ERROR: attrParams expecting 'name = value' pairs separated by comma (%s(%s))\n", paramName, attrname);
            return -1;
          }
          CONFIG_DEBUG_PRINT(("parse_one_attr_param token '%s'\n", token));        
          if (strcmp(token, "=") != 0) {
            fprintf(stderr, "ERROR: attrParams expecting 'name = value' pairs separated by comma (%s(%s))\n", paramName, attrname);
            return 1;
          }
          if (!(token = get_token(infile, tokenbuf))) {
            fprintf(stderr, "ERROR: Did not find value for attrParams %s(%s)\n",
                    paramName, attrname);
            return -1;
          }
          CONFIG_DEBUG_PRINT(("parse_one_attr_param token '%s'\n", token));        
          value = strtod(token, &endptr);
          if (*endptr != '\0') {
            fprintf(stderr, "ERROR: expecting floating point value for attrParam %s(%s) but got '%s'\n", paramName, attrname, token);
            return 1;
          }
          CONFIG_DEBUG_PRINT(("attrParam value %g\n", value));
        }
          
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
        pconfig->attr_param_exponents = (double *)
          safe_realloc(pconfig->attr_param_exponents,
                       (pconfig->num_attr_change_stats_funcs + 1) *
                       sizeof(double));
        pconfig->attr_param_names[pconfig->num_attr_change_stats_funcs] = paramName;
        pconfig->attr_change_stats_funcs[pconfig->num_attr_change_stats_funcs] =
          attr_change_stats_func;
        pconfig->attr_names[pconfig->num_attr_change_stats_funcs] = safe_strdup(attrname);
        pconfig->attr_param_exponents[pconfig->num_attr_change_stats_funcs] =
          exponent_value;
        if (requireErgmValue) {
          pconfig->attr_param_values =
            (double *)safe_realloc(pconfig->attr_param_values,
                  (pconfig->num_attr_change_stats_funcs + 1) * sizeof(double));
          pconfig->attr_param_values[pconfig->num_attr_change_stats_funcs]
            = value;
        }
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
 *
 * If requireErgmValue is TRUE then ERGM parameters require a value
 * (supplied with name = value format) for use in simulation (otherwise
 * no value allowed, used for estimation).
 *
 */
static int parse_attr_params(FILE *infile, param_config_t *pconfig,
                             bool requireErgmValue)
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
                               infile, pconfig, requireErgmValue,
                               ATTR_PARAMS[i].attr_exponent_type)) {
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
 *
 * If requireErgmValue is TRUE then ERGM parameters require a value
 * for use in simulation (otherwise no value allowed, used for
 * estimation). Note that this must be in a different format from that
 * used for normal attribute paramters, as the dyadic parameter as a list
 * (of two) attributes inside the parentheses, but just a single value, e.g.
 * we have 
 *         dyadicParams =  {logGeoDistance(lat,longi) = -0.9513526}
 * rather than for example 
 *         attrParams = {Receiver(teaching_hospital = 0.3689427)}
 *
 * Return nonzero on error else zero.
 */
static int parse_one_dyadic_param(const char *paramName,
                                dyadic_change_stats_func_t *dyadic_change_stats_func,
                                  FILE *infile, param_config_t *pconfig,
                                  bool requireErgmValue)
{
  char      tokenbuf[TOKSIZE];
  char     *token;
  bool      last_token_was_attrname = FALSE;
  bool      opening = TRUE; /* true for first iteration only to expect '(' */
  double    value   = 0;
  char     *endptr; /* for strtod() */
  uint_t    old_num_dyadic = pconfig->num_dyadic_change_stats_funcs;
  uint_t    i;
  
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
  
  if (requireErgmValue) {
    if (!(token = get_token(infile, tokenbuf))) {
      fprintf(stderr, "ERROR: dyadicParams expecting 'dyadicparamname(attr,attr) = value' separated by comma (%s)\n", paramName);
      return -1;
    }
    CONFIG_DEBUG_PRINT(("parse_one_attr_param token '%s'\n", token));        
    if (strcmp(token, "=") != 0) {
      fprintf(stderr, "ERROR: dyadicParams expecting 'dyadicparamname(attr,attr) = value' separated by comma (%s)\n", paramName);
      return 1;
    }
    if (!(token = get_token(infile, tokenbuf))) {
      fprintf(stderr, "ERROR: Did not find value for dyadicParams %s)\n",
              paramName);
      return -1;
    }
    CONFIG_DEBUG_PRINT(("parse_one_dyadic_param token '%s'\n", token));        
    value = strtod(token, &endptr);
    if (*endptr != '\0') {
      fprintf(stderr, "ERROR: expecting floating point value for dyadicParam %s but got '%s'\n", paramName, token);
      return 1;
    }
    CONFIG_DEBUG_PRINT(("dyadicParam value %g\n", value));

    /* Because is has multiple attribute names, we get multiple entries for the
       change statistic - so put the a copy of the value for each one.
       This will be fixed up later in build_dyadic_indices_from_names(). */
    /* TODO this system of ad hoc handling of dyadic parameters was
       alright when there was only one (or now two) types but it will
       be a real mess for any more, should make it more general (and
       hopefully simpler) */
    for (i = old_num_dyadic; i < pconfig->num_dyadic_change_stats_funcs; i++) {
      pconfig->dyadic_param_values =
        (double *)safe_realloc(pconfig->dyadic_param_values,
                               (i + 1) * sizeof(double));
      pconfig->dyadic_param_values[i] = value;
    }
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
 * If requireErgmValue is TRUE then ERGM parameters require a value
 * (supplied with name = value format) for use in simulation (otherwise
 * no value allowed, used for estimation).
 *
 */
static int parse_dyadic_params(FILE *infile, param_config_t *pconfig,
                               bool requireErgmValue)
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
                                 infile, pconfig, requireErgmValue)) {
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
 * If requireErgmValue is TRUE then ERGM parameters require a value
 * (supplied with name = value format) for use in simulation (otherwise
 * no value allowed, used for estimation).
 *
 * For estimation, these parameters have a syntax like this example:
 *
 *   attrInteractionParams = { BinaryPairInteraction(gender_Female, industry_Personal_Goods),
 *                             BinaryPairInteraction(gender_Male, sector_Oil_and_Gas)
 *                           }
 *
 * and for simulation, like this example:
 *
 *  attrInteractionParams = { BinaryPairInteraction(gender_Female, industry_Personal_Goods = 0.1042966),
 *                            BinaryPairInteraction(gender_Male, sector_Oil_and_Gas = 0.4022389)
 *                          }
 *
 * Note that the value applies to the whole parameter (including pair
 * of attributes), despite the syntax making it look like it belongs
 * to the second attribute - this is just an artifact of the way the
 * parser was retrofitted to include values.

 */
static int parse_one_attr_interaction_param(const char *paramName,
                  attr_interaction_change_stats_func_t *attr_interaction_change_stats_func,
                                            FILE *infile, param_config_t *pconfig,
                                            bool requireErgmValue)
{
  char      tokenbuf[TOKSIZE];
  char     *token;
  bool      last_token_was_attrname = FALSE;
  bool      opening = TRUE; /* true for first iteration only to expect '(' */
  uint_t    num_attr_names = 0;
  double    value = 0;
  char     *endptr; /* for strtod() */

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

          if (requireErgmValue) {
            if (!(token = get_token(infile, tokenbuf))) {
              fprintf(stderr, "ERROR: attrInteractionParams expecting 'name1, name2 = value' pairs separated by comma (%s(%s, ?)\n", paramName, pconfig->attr_interaction_pair_names[pconfig->num_attr_interaction_change_stats_funcs].first);
              return -1;
            }
            CONFIG_DEBUG_PRINT(("parse_one_attr_interaction_param token '%s'\n", token));        
            if (strcmp(token, "=") != 0) {
              fprintf(stderr, "ERROR: attrInteractionParams expecting 'name1, name2 = value' pairs separated by comma (%s(%s, ?)\n", paramName, pconfig->attr_interaction_pair_names[pconfig->num_attr_interaction_change_stats_funcs].first);
              return 1;
            }
            if (!(token = get_token(infile, tokenbuf))) {
              fprintf(stderr, "ERROR: Did not find value for attrInteractionParams %s\n",
                      paramName);
              return -1;
            }
            CONFIG_DEBUG_PRINT(("parse_one_attr_interaction_param token '%s'\n", token));        
            value = strtod(token, &endptr);
            if (*endptr != '\0') {
              fprintf(stderr, "ERROR: expecting floating point value for attrInteractionParams %s(%s, %s) but got '%s'\n", paramName, pconfig->attr_interaction_pair_names[pconfig->num_attr_interaction_change_stats_funcs].first, pconfig->attr_interaction_pair_names[pconfig->num_attr_interaction_change_stats_funcs].second, token);
              return 1;
            }
          CONFIG_DEBUG_PRINT(("attrInteractionParam value %g\n", value));

          pconfig->attr_interaction_param_values =
            (double *)safe_realloc(pconfig->attr_interaction_param_values,
                                   (pconfig->num_attr_interaction_change_stats_funcs + 1) * sizeof(double));
          pconfig->attr_interaction_param_values[pconfig->num_attr_interaction_change_stats_funcs]
            = value;
          }
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
 *
 * If requireErgmValue is TRUE then ERGM parameters require a value
 * (supplied with name = value format) for use in simulation (otherwise
 * no value allowed, used for estimation).
 */
static int parse_attr_interaction_params(FILE *infile, param_config_t *pconfig,
                                         bool requireErgmValue)
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
                                           infile, pconfig, requireErgmValue)) {
        fprintf(stderr, "ERROR parsing attrInteractionParams %s\n",
                ATTR_INTERACTION_PARAMS[i].name);
        return 1;
      }
          }
    token = get_token(infile, tokenbuf);
  }
  return 0;
}




/*****************************************************************************
 *
 * external functions
 *
 ****************************************************************************/


/* 
 * Skip whitespace in (open read) file f. First non-whitespace char
 * is returned (or EOF on end or error).
 */
int fskip(FILE *f)
{
  int c = fgetc(f);
  while (c != EOF && isspace(c))
    c = fgetc(f);
  return c;
}

/* Return nonzero if c is a standalone token character else 0 */
int isSingleCharToken(int c)
{
  return c == '=' || c == ',' ||
    c == OPEN_PAREN_CHAR || c == CLOSE_PAREN_CHAR ||
    c == OPEN_SET_CHAR || c == CLOSE_SET_CHAR;
}

/* Return nonzero if c is a character allowed in a token,
   but not a token on its own, else 0 */
int istokenchar(int c)
{
  return c != COMMENT_CHAR && (isalnum(c) || ispunct(c)) &&
    !isSingleCharToken(c);
}

/* return nonzero if c is a char allowed in a parameter name else 0 */
int isParamNameChar(int c)
{
  return isalnum(c) || c == '_';
}

/* return nonzero if s is a lexically valid config parameter name else 0 */
int isValidParamname(const char *s)
{
  const char *p;
  for (p = s; *p != '\0'; p++)
    if (!isParamNameChar(*p))
      return 0;
  return 1;
}

/*
 * get next token in infile (open for read). Returns NULL if none. The
 * token is put in the token parameter (which is return value). token
 * must be allocated by caller with space for at least TOKSIZE chars.
 */
char *get_token(FILE *infile, char *token)
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
  if (c == QUOTE_STR_CHAR) { /* quoted string e.g. "string literal" */
    c = fgetc(infile);
    while (c != EOF && c != QUOTE_STR_CHAR && i < (int)TOKSIZE-1) {
      token[i++] = c;
      token[i] = '\0';
      c = fgetc(infile);
    }
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
int get_paramname_value(FILE *infile, char *paramname, char *value)
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
 * Given name of attribute parameter, returns its type, by searching
 * in the ATTR_PARAMS const array. Note this is a linear search but
 * does not matter as a small hand-coded array anyway.
 */
attr_type_e get_attr_param_type(const char *attrParamName)
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
dyadic_type_e get_dyadic_param_type(const char *dyadicParamName)
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
attr_type_e get_attr_interaction_param_type(const char
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
 * Write the allowed ERGM estimation parameters to stderr 
 */
void dump_parameter_names(void)
{
  uint_t i;
  char   default_lambda_str[128];

  snprintf(default_lambda_str, sizeof(default_lambda_str), " [allows lambda, default %.2f]", DEFAULT_LAMBDA);

  fprintf(stderr, "Structural parameters (%s):\n", STRUCT_PARAMS_STR);
  for (i = 0; i < NUM_STRUCT_PARAMS; i++) {
    fprintf(stderr, "  %s (%s, %s)%s\n", STRUCT_PARAMS[i].name,
            network_type_str(STRUCT_PARAMS[i].network_type),
            network_mode_str(STRUCT_PARAMS[i].network_mode),
            STRUCT_PARAMS[i].struct_param_type == STRUCT_PARAM_TYPE_LAMBDA ?
              default_lambda_str : "");
  }
  fprintf(stderr, "Attribute parameters (%s):\n", ATTR_PARAMS_STR);
  for (i = 0; i < NUM_ATTR_PARAMS; i++) {
    fprintf(stderr, "  %s (%s, %s) (%s)%s\n", ATTR_PARAMS[i].name,
            network_type_str(ATTR_PARAMS[i].network_type),
            network_mode_str(ATTR_PARAMS[i].network_mode),
            ATTR_PARAMS[i].type == ATTR_TYPE_BINARY ? "binary" :
            (ATTR_PARAMS[i].type == ATTR_TYPE_CATEGORICAL ? "categorical" :
             (ATTR_PARAMS[i].type == ATTR_TYPE_CONTINUOUS ? "continuous" :
              (ATTR_PARAMS[i].type == ATTR_TYPE_SET ? "set" :
               "*UNKNOWN*"))),
            ATTR_PARAMS[i].attr_exponent_type == ATTR_EXP_TYPE_EXPONENT ?
            " requires exponent value in [0, 1]" : "");
  }
  fprintf(stderr, "Dyadic covariate parameters (%s):\n", DYADIC_PARAMS_STR);
  for (i = 0; i < NUM_DYADIC_PARAMS; i++) {
    fprintf(stderr, " %s (%s, %s) (%s)\n", DYADIC_PARAMS[i].name,
            network_type_str(DYADIC_PARAMS[i].network_type),
            network_mode_str(DYADIC_PARAMS[i].network_mode),
            DYADIC_PARAMS[i].type == DYADIC_TYPE_GEODISTANCE ?
            "latitude,longitude" :
            (DYADIC_PARAMS[i].type == DYADIC_TYPE_EUCLIDEANDISTANCE ?
             "x, y, z" : "*UNKNOWN*"));
  }
  fprintf(stderr, "Attribute interaction parameters (%s)\n",
          ATTR_INTERACTION_PARAMS_STR);
  for (i = 0; i < NUM_ATTR_INTERACTION_PARAMS; i++) {
    fprintf(stderr, " %s (%s,%s) (%s)\n", ATTR_INTERACTION_PARAMS[i].name,
            network_type_str(ATTR_INTERACTION_PARAMS[i].network_type),
            network_mode_str(ATTR_INTERACTION_PARAMS[i].network_mode),
            ATTR_INTERACTION_PARAMS[i].type == ATTR_TYPE_BINARY ? "binary" :
            (ATTR_INTERACTION_PARAMS[i].type == ATTR_TYPE_CATEGORICAL ? "categorical" :
             (ATTR_INTERACTION_PARAMS[i].type == ATTR_TYPE_CONTINUOUS ? "continuous" :
              "*UNKNOWN*")));
  }
  fprintf(stderr, "\n");
}




/*
 * build_attr_indices_from_names() is called after the config file
 * is parsed by parse_config_file() and also after the 
 * attributes are loaded by calling load_attributes()
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
int build_attr_indices_from_names(param_config_t *pconfig, const graph_t *g)
{
  uint_t i, j;
  bool   found;

  pconfig->attr_indices = safe_malloc(pconfig->num_attr_change_stats_funcs *
                                     sizeof(uint_t));
  
  for (i = 0; i < pconfig->num_attr_change_stats_funcs; i++) {
    found = FALSE;
    switch (get_attr_param_type(pconfig->attr_param_names[i])) {
      case ATTR_TYPE_BINARY:
        for (j = 0; j < g->num_binattr; j++) {
          if (strcasecmp(pconfig->attr_names[i], g->binattr_names[j]) == 0) {
            found = TRUE;
            pconfig->attr_indices[i] = j;
            CONFIG_DEBUG_PRINT(("binattr %s(%s) index %u\n",
                                pconfig->attr_param_names[i],
                                pconfig->attr_names[i], j));
          }
        }
        if (!found) {
          fprintf(stderr, "ERROR: binary attribute %s not found\n",
                  pconfig->attr_names[i]);
          return 1;
        }
        break;
        
      case ATTR_TYPE_CATEGORICAL:
        for (j = 0; j < g->num_catattr; j++) {
          if (strcasecmp(pconfig->attr_names[i], g->catattr_names[j]) == 0) {
            found = TRUE;
            pconfig->attr_indices[i] = j;
            CONFIG_DEBUG_PRINT(("catattr %s(%s) index %u\n",
                                pconfig->attr_param_names[i],
                                pconfig->attr_names[i], j));
          }
        }
        if (!found) {
          fprintf(stderr, "ERROR: categorical attribute %s not found\n",
                  pconfig->attr_names[i]);
          return 1;
        }
        break;

      case ATTR_TYPE_CONTINUOUS:
        for (j = 0; j < g->num_contattr; j++) {
          if (strcasecmp(pconfig->attr_names[i], g->contattr_names[j]) == 0) {
            found = TRUE;
            pconfig->attr_indices[i] = j;
            CONFIG_DEBUG_PRINT(("contattr %s(%s) index %u\n",
                                pconfig->attr_param_names[i],
                                pconfig->attr_names[i], j));
          }
        }
        if (!found) {
          fprintf(stderr, "ERROR: continuous attribute %s not found\n",
                  pconfig->attr_names[i]);
          return 1;
        }
        break;

      case ATTR_TYPE_SET:
        for (j = 0; j < g->num_setattr; j++) {
          if (strcasecmp(pconfig->attr_names[i], g->setattr_names[j]) == 0) {
            found = TRUE;
            pconfig->attr_indices[i] = j;
            CONFIG_DEBUG_PRINT(("setattr %s(%s) index %u\n",
                                pconfig->attr_param_names[i],
                                pconfig->attr_names[i], j));
          }
        }
        if (!found) {
          fprintf(stderr, "ERROR: set attribute %s not found\n",
                  pconfig->attr_names[i]);
          return 1;
        }
        break;
        
      default:
        fprintf(stderr, "ERROR (internal): unknown attribute type %u\n",
                get_attr_param_type(pconfig->attr_param_names[i]));
        return 1;
        break;
    }
  }
  return 0;
}

/*
 * build_dyadic_indices_from_names() is called after the config file
 * is parsed by parse_config_file() and also after the 
 * attributes are loaded by calling load_attributes()
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
int build_dyadic_indices_from_names(param_config_t *pconfig,  graph_t *g,
                                    bool requireErgmValue)
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
  double tmpGeoValue = 0, tmpEuclideanValue = 0;
  
  pconfig->dyadic_indices = safe_malloc(pconfig->num_dyadic_change_stats_funcs *
                                     sizeof(uint_t));
  pconfig->dyadic_types = safe_malloc(pconfig->num_dyadic_change_stats_funcs *
                                     sizeof(uint_t));  
  
  for (i = 0; i < pconfig->num_dyadic_change_stats_funcs; i++) {
    found = FALSE;
    dyadicType = get_dyadic_param_type(pconfig->dyadic_param_names[i]);
    switch (dyadicType) {

      case DYADIC_TYPE_GEODISTANCE:
      case DYADIC_TYPE_EUCLIDEANDISTANCE:

        for (j = 0; j < g->num_contattr; j++) {
          if (strcasecmp(pconfig->dyadic_names[i], g->contattr_names[j]) == 0) {
            found = TRUE;
            pconfig->dyadic_indices[i] = j;
            pconfig->dyadic_types[i] = dyadicType;
            CONFIG_DEBUG_PRINT(("dyadic covariate type %s "
                                "contattr %s(%s) index %u\n",
                                (dyadicType == DYADIC_TYPE_GEODISTANCE ?
                                 "GEODISTANCE" :
                                 (dyadicType == DYADIC_TYPE_EUCLIDEANDISTANCE ?
                                  "EUCLIDEANDISTANCE" : "*ERROR*")),
                                pconfig->dyadic_param_names[i],
                                pconfig->dyadic_names[i], j));
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
                  pconfig->dyadic_names[i]);
          return 1;
        }
        break;

      default:
        fprintf(stderr, "ERROR (internal): unknown dyadic covariate type %u\n",
                get_dyadic_param_type(pconfig->dyadic_param_names[i]));
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
      switch (pconfig->dyadic_types[j]) {
        case DYADIC_TYPE_GEODISTANCE:
          switch (geoIndex) {
            case 0:
              g->latitude_index = pconfig->dyadic_indices[j];
              break;
            case 1:
              g->longitude_index = pconfig->dyadic_indices[j];
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
              g->x_index = pconfig->dyadic_indices[j];
              break;
            case 1:
              g->y_index = pconfig->dyadic_indices[j];
              break;
            case 2:
              g->z_index = pconfig->dyadic_indices[j];
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
                  (int)pconfig->dyadic_types[j]);
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
      assert(get_dyadic_param_type(pconfig->dyadic_param_names[firstGeoIndex]) == DYADIC_TYPE_GEODISTANCE);
      assert(get_dyadic_param_type(pconfig->dyadic_param_names[firstEuclideanIndex]) == DYADIC_TYPE_EUCLIDEANDISTANCE);
      assert(pconfig->num_dyadic_change_stats_funcs == 5);
      tmpGeoName = safe_strdup(pconfig->dyadic_names[firstGeoIndex]);
      tmpGeoParamName = safe_strdup(pconfig->dyadic_param_names[firstGeoIndex]);
      tmpGeoFunc = pconfig->dyadic_change_stats_funcs[firstGeoIndex];
      tmpEuclideanName = safe_strdup(pconfig->dyadic_names[firstEuclideanIndex]);
      tmpEuclideanParamName = safe_strdup(pconfig->dyadic_param_names[firstEuclideanIndex]);
      tmpEuclideanFunc = pconfig->dyadic_change_stats_funcs[firstEuclideanIndex];
      if (requireErgmValue) {
        tmpGeoValue = pconfig->dyadic_param_values[firstGeoIndex];
        tmpEuclideanValue = pconfig->dyadic_param_values[firstEuclideanIndex];
      }
      pconfig->dyadic_names[0] = tmpGeoName;
      pconfig->dyadic_param_names[0] = tmpGeoParamName;
      pconfig->dyadic_change_stats_funcs[0] = tmpGeoFunc;
      pconfig->dyadic_types[0] = DYADIC_TYPE_GEODISTANCE;
      pconfig->dyadic_names[1] = tmpEuclideanName;
      pconfig->dyadic_param_names[1] = tmpEuclideanParamName;
      pconfig->dyadic_change_stats_funcs[1] = tmpEuclideanFunc;
      pconfig->dyadic_types[1] = DYADIC_TYPE_EUCLIDEANDISTANCE;
      pconfig->num_dyadic_change_stats_funcs = 2;      
      free(pconfig->dyadic_names[2]);
      free(pconfig->dyadic_names[3]);
      free(pconfig->dyadic_names[4]);
      if (requireErgmValue) {
        pconfig->dyadic_param_values[0] = tmpGeoValue;
        pconfig->dyadic_param_values[1] = tmpEuclideanValue;
      }
    } else if (numGeoAttr > 0) {
      /* only [log]GeoDistance */
      assert(firstGeoIndex >= 0);
      assert(pconfig->num_dyadic_change_stats_funcs == 2);
      assert(get_dyadic_param_type(pconfig->dyadic_param_names[0])
             == DYADIC_TYPE_GEODISTANCE &&
             get_dyadic_param_type(pconfig->dyadic_param_names[1])
             == DYADIC_TYPE_GEODISTANCE);
      pconfig->num_dyadic_change_stats_funcs = 1;
      free(pconfig->dyadic_names[1]);
    } else if (numEuclideanAttr > 0) {
      /* only EuclideanDistance */
      assert(firstEuclideanIndex >= 0);
      assert(pconfig->num_dyadic_change_stats_funcs == 3);
      assert(get_dyadic_param_type(pconfig->dyadic_param_names[0])
             == DYADIC_TYPE_EUCLIDEANDISTANCE &&
             get_dyadic_param_type(pconfig->dyadic_param_names[1])
             == DYADIC_TYPE_EUCLIDEANDISTANCE &&
             get_dyadic_param_type(pconfig->dyadic_param_names[2])
             == DYADIC_TYPE_EUCLIDEANDISTANCE);
      pconfig->num_dyadic_change_stats_funcs = 1;
      free(pconfig->dyadic_names[1]);
      free(pconfig->dyadic_names[2]);
    } else {
      assert(FALSE);
    }
    
  }
  return 0;
}


/*
 * build_attr_interaction_pair_indices_from_names() is called after the
 * config file is parsed by parse_config_file() and also after the
 * attributes are loaded by calling load_attributes()
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
int build_attr_interaction_pair_indices_from_names(param_config_t *pconfig,
                                                   const graph_t *g)
{
  uint_t i, j, attrnum;
  bool   found;
  char  *attrname;

  pconfig->attr_interaction_pair_indices = safe_malloc(
    pconfig->num_attr_interaction_change_stats_funcs * sizeof(uint_pair_t));

  for (i = 0; i < pconfig->num_attr_interaction_change_stats_funcs; i++) {
    for (attrnum = 0; attrnum < 2 /* 0=first and 1=second */; attrnum++) {
      found = FALSE;
      if (attrnum == 0) {
        attrname = pconfig->attr_interaction_pair_names[i].first;
      } else {
        attrname = pconfig->attr_interaction_pair_names[i].second;
      }
      switch (get_attr_interaction_param_type(
                pconfig->attr_interaction_param_names[i])) {
        case ATTR_TYPE_BINARY:
          for (j = 0; j < g->num_binattr; j++) {
            if (strcasecmp(attrname, g->binattr_names[j]) == 0) {
              found = TRUE;
              if (attrnum == 0) {
                pconfig->attr_interaction_pair_indices[i].first = j;
              } else {
                pconfig->attr_interaction_pair_indices[i].second = j;
              }
              CONFIG_DEBUG_PRINT(("binattr interaction [%s] %s(%s) index %u\n",
                                  attrnum == 0 ? "first" : "second",
                                  pconfig->attr_interaction_param_names[i],
                                  attrname, j));
            }
          }
          if (!found) {
            fprintf(stderr, "ERROR: binary attribute %s not found\n",
                    attrname);
            return 1;
          }
          break;

        case ATTR_TYPE_CATEGORICAL:
          for (j = 0; j < g->num_catattr; j++) {
            if (strcasecmp(attrname, g->catattr_names[j]) == 0) {
              found = TRUE;
              if (attrnum == 0) {
                pconfig->attr_interaction_pair_indices[i].first = j;
              } else {
                pconfig->attr_interaction_pair_indices[i].second = j;
              }
              CONFIG_DEBUG_PRINT(("catattr interaction [%s] %s(%s) index %u\n",
                                  attrnum == 0 ? "first" : "second",
                                  pconfig->attr_interaction_param_names[i],
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
 * Free the param config structure
 *
 * Parameters:
 *     pconfig - pointer to param config struct
 *
 * Return value:
 *     None
 */
void free_param_config_struct(param_config_t *pconfig)
{
  uint_t i;

  free(pconfig->change_stats_funcs);
  free(pconfig->param_names);
  free(pconfig->param_values);
  free(pconfig->attr_param_names);
  for (i = 0; i < pconfig->num_attr_change_stats_funcs; i++) 
    free(pconfig->attr_names[i]);
  free(pconfig->attr_indices);
  for (i = 0; i < pconfig->num_dyadic_change_stats_funcs; i++) 
    free(pconfig->dyadic_names[i]);
  free(pconfig->dyadic_indices);
  free(pconfig->dyadic_types);
  free(pconfig->dyadic_param_values);
  free(pconfig->attr_interaction_param_names);
  for (i = 0; i < pconfig->num_attr_interaction_change_stats_funcs; i++)  {
    free(pconfig->attr_interaction_pair_names[i].first);
    free(pconfig->attr_interaction_pair_names[i].second);
  }
  free(pconfig->attr_interaction_pair_indices);
}



/*
 * Write the allowed configuration parameters, their descriptions and 
 * default values, to stderr
 */
void dump_config_names(const void *config,
                       const config_param_t *config_params,
                       uint_t num_config_params)
{
  uint_t i;
  fprintf(stderr, "Configuration parameters:\n");
  for (i = 0; i < num_config_params; i++) {
    fprintf(stderr, "  %s: %s ", config_params[i].name,
            config_params[i].description);
    switch (config_params[i].type) {
      case PARAM_TYPE_DOUBLE:
        fprintf(stderr, "(floating point) [default %g]\n",
                *(const double *)((const char *)config + config_params[i].offset));
        break;
        
      case PARAM_TYPE_UINT:
        fprintf(stderr, "(unsigned integer) [default %u]\n",
                *(const uint_t *)((const char *)config + config_params[i].offset));
        break;

      case PARAM_TYPE_BOOL:
        fprintf(stderr, "(Boolean) [default %s]\n",
                *(const bool *)((const char *)config + config_params[i].offset) ?
                "True" : "False");
        break;

      case PARAM_TYPE_STRING:
        fprintf(stderr, "(string)");
        if (*(const char * const *)((const char *)config + config_params[i].offset))
          fprintf(stderr, " [default %s]\n", *(const char * const *)((const char *)config + config_params[i].offset));
        else
          fprintf(stderr, "\n");
        break;

      case PARAM_TYPE_SET:
        fprintf(stderr, "(set of ERGM parameter names)\n");
        break;

      case PARAM_TYPE_ULONG:
        fprintf(stderr, "(unsigned long integer) [default %lu]\n",
                *(const ulong_t *)((const char *)config + config_params[i].offset));
        break;

        
      default:
      fprintf(stderr, "ERROR (internal): unknown parameter type %d\n",
              config_params[i].type);
      break;
    }
  }
}


/* 
 * Given the parameter name and value parsed from the config file,
 * check that they are valid, using the file config_params (constant
 * input param), and set the corresponding fields in the config and
 * pconfig (output param) structures, using hte config_is_set bool
 * array (set to all FALSE initially, in/out parameter) to check for
 * duplicates. Returns nonzero on error, else zero.  The infile (open
 * read) from which they are parsed is also passed.  For set values,
 * the valuestr is just '{' and more parsing is required, for which
 * the infile is passed to the appropriate function.
 * If requireErgmValue is TRUE then ERGM parameters require a value
 * (supplied with name = value format) for use in simulation (otherwise
 * no value allowed, used for estimation).
 */
int check_and_set_param_value(const char *paramname,
                              const char *valuestr,
                              FILE *infile,
                              void *config,
                              bool *config_is_set,
                              param_config_t *pconfig,
                              const config_param_t *config_params,
                              uint_t num_config_params,
                              bool requireErgmValue)
{
  uint_t i;
  bool   found_paramname = FALSE;
  double valuedouble;
  char  *endptr; /* for strtod() */
  uint_t valueint;
  bool   valuebool;
  ulong_t valuelong;

  for (i = 0; i < num_config_params; i++) {
    if (strcasecmp(paramname, config_params[i].name) == 0) {
      found_paramname = TRUE;
      break;
    }
  }
  if (!found_paramname) {
    fprintf(stderr, "ERROR: invalid parameter name '%s'\n", paramname);
    return 1;
  }

  if (config_is_set[i]) {
    fprintf(stderr, "ERROR: parameter %s is set more than once\n",
            config_params[i].name);
    return 1;
  }
  config_is_set[i] = TRUE;
  
  switch(config_params[i].type) {
    case PARAM_TYPE_DOUBLE:  /* numeric (floating point) */
      valuedouble = strtod(valuestr, &endptr);
      if (*endptr != '\0') {
        fprintf(stderr, "ERROR: expecting floating point value for parameter %s but got '%s'\n", paramname, valuestr);
        return 1;
      }
      CONFIG_DEBUG_PRINT(("%s = %g\n", paramname, valuedouble));
      *(double *)((char *)config + config_params[i].offset) = valuedouble;
      break;
      
    case PARAM_TYPE_UINT:    /* numeric (unsigned integer) */
      if (sscanf(valuestr, "%u", &valueint) != 1) {
        fprintf(stderr, "ERROR: expecting unsigned integer value for parameter %s but got '%s'\n", paramname, valuestr);
        return 1;
      }
      CONFIG_DEBUG_PRINT(("%s = %u\n", paramname, valueint));
      *(uint_t *)((char *)config + config_params[i].offset) = valueint;
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
      *(bool *)((char *)config + config_params[i].offset) = valuebool;
      break;
      
    case PARAM_TYPE_STRING:  /* string (may be quoted, not necessarily) */
      /* free default/previous value first (may be NULL, which is alright)*/
      free(*(char **)((char *)config + config_params[i].offset));
      *(char **)((char *)config + config_params[i].offset) =
        safe_strdup(valuestr);
      CONFIG_DEBUG_PRINT(("%s = %s\n", paramname, valuestr));
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
        if (pconfig->num_change_stats_funcs > 0) {
          fprintf(stderr, "ERROR: %s specified more than once\n",
                  STRUCT_PARAMS_STR);
          return 1;
        } 
        return parse_struct_params(infile, pconfig, requireErgmValue); 
      } else if (strcasecmp(paramname, ATTR_PARAMS_STR) == 0) {
        if (pconfig->num_attr_change_stats_funcs > 0) {
          fprintf(stderr, "ERROR: %s specified more than once\n",
                  ATTR_PARAMS_STR);
          return 1;
        }
        return parse_attr_params(infile, pconfig, requireErgmValue);
      } else if (strcasecmp(paramname, DYADIC_PARAMS_STR) == 0) {        
        if (pconfig->num_dyadic_change_stats_funcs > 0) {
          fprintf(stderr, "ERROR: %s specified more than once\n",
                  DYADIC_PARAMS_STR);
          return 1;
        }
        return parse_dyadic_params(infile, pconfig, requireErgmValue);
      } else if (strcasecmp(paramname, ATTR_INTERACTION_PARAMS_STR) == 0) {
        if (pconfig->num_attr_interaction_change_stats_funcs > 0) {
          fprintf(stderr, "ERROR: %s specified more than once\n",
                  ATTR_INTERACTION_PARAMS_STR);
          return 1;
        }
        return parse_attr_interaction_params(infile, pconfig, requireErgmValue);
      } else {
        fprintf(stderr, "ERROR (internal): unknown parameter %s\n", paramname);
        return 1;
      }
      break;

    case PARAM_TYPE_ULONG:    /* numeric (unsigned long integer) */
      if (sscanf(valuestr, "%lu", &valuelong) != 1) {
        fprintf(stderr, "ERROR: expecting unsigned long integer value for parameter %s but got '%s'\n", paramname, valuestr);
        return 1;
      }
      CONFIG_DEBUG_PRINT(("%s = %lu\n", paramname, valuelong));
      *(ulong_t *)((char *)config + config_params[i].offset) = valuelong;
      break;
        
    default:
      fprintf(stderr, "ERROR (internal): unknown parameter type %d\n",
              config_params[i].type);
      return 1;
      break;
  }
  return 0;
}


/*
 * Given name of structural parameter, returns its network type, by searching
 * in the STRUCT_PARAMS const array. Note this is a linear search but
 * does not matter as a small hand-coded array anyway.
 */
network_type_e get_struct_param_network_type(const char *structParamName)
{
  uint_t i;

  for (i = 0; i < NUM_STRUCT_PARAMS; i++) {
    if (strcasecmp(STRUCT_PARAMS[i].name, structParamName) == 0) {
      return STRUCT_PARAMS[i].network_type;
    }
  }
  return NETWORK_TYPE_INVALID;
}



/*
 * Given name of attribute parameter, returns its network type, by searching
 * in the ATTR_PARAMS const array. Note this is a linear search but
 * does not matter as a small hand-coded array anyway.
 */
network_type_e get_attr_param_network_type(const char *attrParamName)
{
  uint_t i;

  for (i = 0; i < NUM_ATTR_PARAMS; i++) {
    if (strcasecmp(ATTR_PARAMS[i].name, attrParamName) == 0) {
      return ATTR_PARAMS[i].network_type;
    }
  }
  return NETWORK_TYPE_INVALID;
}

/*
 * Given name of dyadic covariate parameter, returns its network type, 
 * by searching
 * in the DYADIC_PARAMS const array. Note this is a linear search but
 * does not matter as a small hand-coded array anyway.
 */
network_type_e get_dyadic_param_network_type(const char *dyadicParamName)
{
  uint_t i;

  for (i = 0; i < NUM_DYADIC_PARAMS; i++) {
    if (strcasecmp(DYADIC_PARAMS[i].name, dyadicParamName) == 0) {
      return DYADIC_PARAMS[i].network_type;
    }
  }
  return NETWORK_TYPE_INVALID;
}

/*
 * Given name of attribute interaction parameter, returns its network type,
 * by searching
 * in the ATTR_INTERACTION_PARAMS const array. Note this is a linear search but
 * does not matter as a small hand-coded array anyway.
 */
network_type_e get_attr_interaction_param_network_type(const char
                                                       *attrInteractionParamName)
{
  uint_t i;

  for (i = 0; i < NUM_ATTR_INTERACTION_PARAMS; i++) {
    if (strcasecmp(ATTR_INTERACTION_PARAMS[i].name,
                   attrInteractionParamName) == 0) {
      return ATTR_INTERACTION_PARAMS[i].network_type;
    }
  }
  return NETWORK_TYPE_INVALID;
}



/*
 * Given name of structural parameter, returns its network mode, by searching
 * in the STRUCT_PARAMS const array. Note this is a linear search but
 * does not matter as a small hand-coded array anyway.
 */
network_mode_e get_struct_param_network_mode(const char *structParamName)
{
  uint_t i;

  for (i = 0; i < NUM_STRUCT_PARAMS; i++) {
    if (strcasecmp(STRUCT_PARAMS[i].name, structParamName) == 0) {
      return STRUCT_PARAMS[i].network_mode;
    }
  }
  return NETWORK_MODE_INVALID;
}



/*
 * Given name of attribute parameter, returns its network mode, by searching
 * in the ATTR_PARAMS const array. Note this is a linear search but
 * does not matter as a small hand-coded array anyway.
 */
network_mode_e get_attr_param_network_mode(const char *attrParamName)
{
  uint_t i;

  for (i = 0; i < NUM_ATTR_PARAMS; i++) {
    if (strcasecmp(ATTR_PARAMS[i].name, attrParamName) == 0) {
      return ATTR_PARAMS[i].network_mode;
    }
  }
  return NETWORK_MODE_INVALID;
}

/*
 * Given name of dyadic covariate parameter, returns its network mode, 
 * by searching
 * in the DYADIC_PARAMS const array. Note this is a linear search but
 * does not matter as a small hand-coded array anyway.
 */
network_mode_e get_dyadic_param_network_mode(const char *dyadicParamName)
{
  uint_t i;

  for (i = 0; i < NUM_DYADIC_PARAMS; i++) {
    if (strcasecmp(DYADIC_PARAMS[i].name, dyadicParamName) == 0) {
      return DYADIC_PARAMS[i].network_mode;
    }
  }
  return NETWORK_MODE_INVALID;
}

/*
 * Given name of attribute interaction parameter, returns its network mode,
 * by searching
 * in the ATTR_INTERACTION_PARAMS const array. Note this is a linear search but
 * does not matter as a small hand-coded array anyway.
 */
network_mode_e get_attr_interaction_param_network_mode(const char
                                                       *attrInteractionParamName)
{
  uint_t i;

  for (i = 0; i < NUM_ATTR_INTERACTION_PARAMS; i++) {
    if (strcasecmp(ATTR_INTERACTION_PARAMS[i].name,
                   attrInteractionParamName) == 0) {
      return ATTR_INTERACTION_PARAMS[i].network_mode;
    }
  }
  return NETWORK_MODE_INVALID;
}



/*
 * Return string describing network type
 */
const char *network_type_str(network_type_e net_type)
{
  switch(net_type) {
    case NETWORK_TYPE_DIRECTED:
      return "Directed";
      break;

    case NETWORK_TYPE_UNDIRECTED:
      return "Undirected";
      break;

    case NETWORK_TYPE_BOTH:
      return "Both";

    default:
      return "*UNKNOWN*";
      break;
  }
}

/*
 * Return string describing network mode
 */
const char *network_mode_str(network_mode_e net_mode)
{
  switch(net_mode) {
    case NETWORK_MODE_ONEMODE:
      return "One-mode";
      break;

    case NETWORK_MODE_TWOMODE:
      return "Two-mode";
      break;

    case NETWORK_MODE_BOTH:
      return "Both";

    default:
      return "*UNKNOWN*";
      break;
  }
}

/*
 * Return TRUE if network type (used to label ERGM parameters, i.e. change
 * statistics functions) is allowed with the given graph  else FALSE
 */
bool is_allowed_network_type(network_type_e net_type, const graph_t *g)
{
  switch(net_type) {
    case NETWORK_TYPE_DIRECTED:
      return g->is_directed;
      break;

    case NETWORK_TYPE_UNDIRECTED:
      return !g->is_directed;
      break;

    case NETWORK_TYPE_BOTH:
      return TRUE;

    default:
      fprintf(stderr, "ERROR (internal): Unknown network type %d\n", net_type);
      return FALSE;
      break;
  }
}

/*
 * Return TRUE if network mode (used to label ERGM parameters, i.e. change
 * statistics functions) is allowed with the given graph  else FALSE
 */
bool is_allowed_network_mode(network_mode_e net_mode, const graph_t *g)
{
  switch(net_mode) {
    case NETWORK_MODE_ONEMODE:
      return !g->is_bipartite;
      break;

    case NETWORK_MODE_TWOMODE:
      return g->is_bipartite;
      break;

    case NETWORK_MODE_BOTH:
      return TRUE;

    default:
      fprintf(stderr, "ERROR (internal): Unknown network mode %d\n", net_mode);
      return FALSE;
      break;
  }
}



/*
 * Check all the ERGM paramters are allowed for the network type
 * (directed / undirected) and mode (one-mode / two-mode)
 * Return nonzero on error.
 */
int check_param_network_type(param_config_t *pconfig, const graph_t *g)
{
  uint_t         i;
  network_type_e net_type;
  network_mode_e net_mode;

  /* structural parameters */
  for (i = 0; i < pconfig->num_change_stats_funcs; i++) {
    net_type = get_struct_param_network_type(pconfig->param_names[i]);
    if (!(is_allowed_network_type(net_type, g))) {
      fprintf(stderr, "Bad network type (%s) for %s %s (%s)\n",
              g->is_directed ? "Directed" : "Undirected",
              STRUCT_PARAMS_STR,
              pconfig->param_names[i],
              network_type_str(net_type));

      return 1;
    }
    net_mode = get_struct_param_network_mode(pconfig->param_names[i]);
    if (!(is_allowed_network_mode(net_mode, g))) {
      fprintf(stderr, "Bad network mode (%s) for %s %s (%s)\n",
              g->is_bipartite ? "Two-mode" : "One-mode",
              STRUCT_PARAMS_STR,
              pconfig->param_names[i],
              network_mode_str(net_mode));

      return 1;
    }
  }

  /* attribute parameters */
  for (i = 0; i < pconfig->num_attr_change_stats_funcs; i++) {
    net_type = get_attr_param_network_type(pconfig->attr_param_names[i]);
    if (!(is_allowed_network_type(net_type, g))) {
      fprintf(stderr, "Bad network type (%s) for %s %s (%s)\n",
              g->is_directed ? "Directed" : "Undirected",
              ATTR_PARAMS_STR,
              pconfig->attr_param_names[i],
              network_type_str(net_type));

      return 1;
    }
    net_mode = get_attr_param_network_mode(pconfig->attr_param_names[i]);
    if (!(is_allowed_network_mode(net_mode, g))) {
      fprintf(stderr, "Bad network mode (%s) for %s %s (%s)\n",
              g->is_bipartite ? "Two-mode" : "One-mode",
              ATTR_PARAMS_STR,
              pconfig->attr_param_names[i],
              network_mode_str(net_mode));

      return 1;
    }
  }

  /* dyadic parameters */
  for (i = 0; i < pconfig->num_dyadic_change_stats_funcs; i++) {
    net_type = get_dyadic_param_network_type(pconfig->dyadic_param_names[i]);
    if (!(is_allowed_network_type(net_type, g))) {
      fprintf(stderr, "Bad network type (%s) for %s %s (%s)\n",
              g->is_directed ? "Directed" : "Undirected",
              DYADIC_PARAMS_STR,
              pconfig->dyadic_param_names[i],
              network_type_str(net_type));
      return 1;
    }
    net_mode = get_dyadic_param_network_mode(pconfig->dyadic_param_names[i]);
    if (!(is_allowed_network_mode(net_mode, g))) {
      fprintf(stderr, "Bad network mode (%s) for %s %s (%s)\n",
              g->is_bipartite ? "Two-mode" : "One-mode",
              DYADIC_PARAMS_STR,
              pconfig->dyadic_param_names[i],
              network_mode_str(net_mode));
      return 1;
    }
  }

  /* attribute interaction parameters */
  for (i = 0; i < pconfig->num_attr_interaction_change_stats_funcs; i++) {
    net_type = get_attr_interaction_param_network_type(pconfig->attr_interaction_param_names[i]);
    if (!(is_allowed_network_type(net_type, g))) {
      fprintf(stderr, "Bad network type (%s) for %s %s (%s)\n",
              g->is_directed ? "Directed" : "Undirected",
              ATTR_INTERACTION_PARAMS_STR,
              pconfig->attr_interaction_param_names[i],
              network_type_str(net_type));
      return 1;
    }
    net_mode = get_attr_interaction_param_network_mode(pconfig->attr_interaction_param_names[i]);
    if (!(is_allowed_network_mode(net_mode, g))) {
      fprintf(stderr, "Bad network mode (%s) for %s %s (%s)\n",
              g->is_bipartite ? "Two-mode" : "One-mode",
              ATTR_INTERACTION_PARAMS_STR,
              pconfig->attr_interaction_param_names[i],
              network_mode_str(net_mode));
      return 1;
    }
  }
  return 0;
}
