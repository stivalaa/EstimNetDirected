/*****************************************************************************
 * 
 * File:    EstimNetDirectedMPImain.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * The MPI version simply runs multiple instances of the estimation
 * (each parsing the same config file) in parallel, writing the output
 * to separpate files (common prefix with MPI rank suffix).
 * 
 *
 *   Usage: EstimNetDirected_mpi config_filename
 *
 ****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <unistd.h>
#include <mpi.h>
#include "utils.h"
#include "estimconfigparser.h"
#include "equilibriumExpectation.h"

/*****************************************************************************
 *
 * Constants
 *
 ****************************************************************************/

static const int MPI_RANK_MASTER  = 0; /* MPI master task rank number */

/*****************************************************************************
 *
 * File static variables
 *
 ****************************************************************************/

static char myname[MPI_MAX_PROCESSOR_NAME];   /* MPI name */
static int  mynamelen;                        /* length of myname */
static int  numtasks, rank;                   /* MPI number of tasks, rank */


/*****************************************************************************
 *
 * Main
 *
 ****************************************************************************/

static void usage(const char *progname)
{
  fprintf(stderr, "Usage: %s [-h] config_filename\n"
          "  -h : write parameter names to stderr and exit\n"
          , progname);
  exit(1);
}

int main(int argc, char *argv[])
{
  int        c;
  char      *config_filename = NULL;
  config_t  *config;
  int        rc;

  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    fprintf (stderr, "Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Get_processor_name(myname, &mynamelen);

  init_prng(rank); /* initialize pseudorandom number generator */

  init_config_parser();

  while ((c = getopt(argc, argv, "h")) != -1)  {
    switch (c)   {
      case 'h':
        dump_config_names(&CONFIG);
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

  printf("MPI name %s rank %d of total %d\n",myname,rank,numtasks);
  
  config_filename = argv[optind];
  if (!(config = parse_config_file(config_filename))) {
    if (rank == MPI_RANK_MASTER) {
      fprintf(stderr, "ERROR parsing configuration file %s\n", config_filename);
    }
    rc = 1;
  } else {
    rc = do_estimation(config, rank);
  }
  free_config_struct(config);
  MPI_Finalize();
  exit(rc);
}

