/*****************************************************************************
 * 
 * File:    EstimNetDirectedMain.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 *
 *   Usage: EstimNetDirected config_filename
 *
 ****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <unistd.h>
#include "utils.h"
#include "configparser.h"
#include "equilibriumExpectation.h"


/*****************************************************************************
 *
 * Main
 *
 ****************************************************************************/

static void usage(const char *progname)
{
  fprintf(stderr, "Usage: %s [-h] config_filename\n"
          "  -h : write parameter names to stderr\n"
          , progname);
  exit(1);
}

int main(int argc, char *argv[])
{
  int        c;
  char      *config_filename = NULL;
  config_t  *config;
  int        rc;

  init_prng(0); /* initialize pseudorandom number generator */

  init_config_parser();
  
  while ((c = getopt(argc, argv, "h")) != -1)  {
    switch (c)   {
      case 'h':
        dump_config_names();
        dump_parameter_names();
        break;
      default:
        usage(argv[0]);
        break;
    }
  }

  if (argc - optind != 1)
    usage(argv[0]);

  config_filename = argv[optind];
  if (!(config = parse_config_file(config_filename))) {
    fprintf(stderr, "ERROR parsing configuration file %s\n", config_filename);
    rc = 1;
  } else {
    rc = do_estimation(config, 0);
  }
  free_config_struct(config);
  exit(rc);
}

