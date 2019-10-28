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
 * constant definitions
 *
 ****************************************************************************/

const size_t BUFSIZE = 16384;  /* line buffer size for reading files */
const size_t TOKSIZE = 1024;   /* maximum size of a token */
const char   COMMENT_CHAR = '#'; /* comment character */
const char   OPEN_SET_CHAR = '{';  /* set of parameter vals open */
const char   CLOSE_SET_CHAR = '}'; /* set of parameter vals close */
const char   OPEN_PAREN_CHAR = '(';
const char   CLOSE_PAREN_CHAR = ')';


/* True and False values for Boolean config value type. Not case sensitive */
const char *TRUE_STR = "true";
const char *FALSE_STR = "false";


/*
 * Structural parameters allowed as the names in the set for the 
 * structParams parameter. Names are not case sensitive 
 */
const struct_param_t STRUCT_PARAMS[] =
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
const uint_t NUM_STRUCT_PARAMS = sizeof(STRUCT_PARAMS) /
  sizeof(STRUCT_PARAMS[0]);

/*
 * Attribute parameters allowed as the names in the set for the 
 * attrParams  parameters. Names are not case sensitive.
 */
const attr_param_t ATTR_PARAMS[] =
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
  {"DiffSign",               ATTR_TYPE_CONTINUOUS, changeDiffSign},
  {"DiffDirSR",              ATTR_TYPE_CONTINUOUS, changeDiffDirSR},
  {"DiffDirRS",              ATTR_TYPE_CONTINUOUS, changeDiffDirRS},
  {"JaccardSimilarity",      ATTR_TYPE_SET,        changeJaccardSimilarity}
};
const uint_t NUM_ATTR_PARAMS = sizeof(ATTR_PARAMS) /
  sizeof(ATTR_PARAMS[0]);


/*
 * Dyadic covariate parameters allowed as the names in the set for the
 * dyadicParams parameters. Names are not case sensitive.
 */
const dyadic_param_t DYADIC_PARAMS[] =
{
  {"GeoDistance",    DYADIC_TYPE_GEODISTANCE,   changeGeoDistance},
  {"logGeoDistance", DYADIC_TYPE_GEODISTANCE,   changeLogGeoDistance},
  {"EuclideanDistance", DYADIC_TYPE_EUCLIDEANDISTANCE, changeEuclideanDistance}
};
const uint_t NUM_DYADIC_PARAMS = sizeof(DYADIC_PARAMS) /
  sizeof(DYADIC_PARAMS[0]);

/*
 * Attribute pair interaction parameters allowed as the names in the
 * set for the attrInteractionParams parameters. Names are not case sensitive.
 */
const attr_interaction_param_t ATTR_INTERACTION_PARAMS[] =
{
  {"MatchingInteraction",     ATTR_TYPE_CATEGORICAL, changeMatchingInteraction},
};
const uint_t NUM_ATTR_INTERACTION_PARAMS =
  sizeof(ATTR_INTERACTION_PARAMS) / sizeof(ATTR_INTERACTION_PARAMS[0]);



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
              (ATTR_PARAMS[i].type == ATTR_TYPE_SET ? "set" :
               "*UNKNOWN*"))));
  }
  fprintf(stderr, "Dyadic covariate parameters (%s):\n", DYADIC_PARAMS_STR);
  for (i = 0; i < NUM_DYADIC_PARAMS; i++) {
    fprintf(stderr, " %s (%s)\n", DYADIC_PARAMS[i].name,
            DYADIC_PARAMS[i].type == DYADIC_TYPE_GEODISTANCE ?
            "latitude,longitude" :
            (DYADIC_PARAMS[i].type == DYADIC_TYPE_EUCLIDEANDISTANCE ?
             "x, y, z" : "*UNKNOWN*"));
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

