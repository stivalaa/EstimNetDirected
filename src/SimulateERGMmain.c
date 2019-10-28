/*****************************************************************************
 * 
 * File:    SimulateERGMmain.c
 * Author:  Alex Stivala
 * Created: October 2019
 *
 *
 *   Usage: SimulateERGM sim_config_filename
 *
 ****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <unistd.h>
#include "utils.h"
#include "simconfigparser.h"


/*****************************************************************************
 *
 * Main
 *
 ****************************************************************************/

static void usage(const char *progname)
{
  fprintf(stderr, "Usage: %s [-h] sim_config_filename\n"
          "  -h : write parameter names to stderr and exit\n"
          , progname);
  exit(1);
}

int main(int argc, char *argv[])
{
  int              c;
  char            *config_filename = NULL;
  sim_config_t    *config;
  int              rc;

  init_prng(0); /* initialize pseudorandom number generator */

  init_sim_config_parser();
  
  while ((c = getopt(argc, argv, "h")) != -1)  {
    switch (c)   {
      case 'h':
        dump_config_names(&SIM_CONFIG, (const config_param_t *)&SIM_CONFIG_PARAMS, NUM_SIM_CONFIG_PARAMS);
        dump_parameter_names();
        exit(0);
        break;
      default:
        usage(argv[0]);
        break;
    }
  }

  if (argc - optind != 1)
    usage(argv[0]);

  config_filename = argv[optind];
  if (!(config = parse_sim_config_file(config_filename))) {
    fprintf(stderr, "ERROR parsing configuration file %s\n", config_filename);
    rc = 1;
  } else {
    rc = do_simulation(config, 0);
  }
  free_sim_config_struct(config);
  exit(rc);
}

